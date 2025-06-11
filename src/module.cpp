#include "module.hh"

#include <iostream>
#include <math.h>
#include <array>
#include <limits>
#include <unordered_set>
#include <unordered_map>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <TinyAD/Utils/Out.hh>

namespace std
{
    template <>
    struct hash<MakePlanarFacesPlus::MeshEdge> {
        size_t operator()(const MakePlanarFacesPlus::MeshEdge& e) const noexcept {
            size_t hash = 0;
            MakePlanarFacesPlus::hash_combine(hash, e.v1);
            MakePlanarFacesPlus::hash_combine(hash, e.v2);
            return hash;
        }
    };
}

namespace py = pybind11;

namespace MakePlanarFacesPlus
{
    inline void hash_combine(std::size_t& hash, const int& v)
    {
        // https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
        std::hash<int> hasher;
        hash ^= hasher(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }

    std::vector<Vec3d> make_planar_faces(const std::vector<Vec3d>& vertices, const std::vector<std::vector<int>>& faces, const std::vector<int>& fixed_vertices, const MakePlanarSettings& settings)
    {
        // Store fixed vertices in hashmap for faster lookup
        std::unordered_set<int> fixed_vertices_set(fixed_vertices.begin(), fixed_vertices.end());

        // Go over all faces and store connected edge segments for planarity objective
        // Also compute the min edge length of the mesh and store all edges
        std::vector<std::array<int, 4>> planarity_objective_indices;
        double min_edge_length = std::numeric_limits<double>::infinity();
        std::unordered_set<MeshEdge> mesh_edge_set;
        int skipped_faces = 0;
        for (std::vector<int> face : faces)
        {
            int n_verts = face.size();
            TINYAD_ASSERT_GEQ(n_verts, 3);
            if (n_verts == 3) continue;
            bool some_vertex_free = false;
            for (const int& v_id : face)
            {
                if (!some_vertex_free && fixed_vertices_set.find(v_id) == fixed_vertices_set.end()) some_vertex_free = true;
            }
            if (!some_vertex_free)
            {
                // Skip a face that is completely pinned in place
                skipped_faces++;
                continue;
            }
            for (int v_id = 0; v_id < n_verts; v_id++)
            {
                planarity_objective_indices.push_back({ face[v_id], face[(v_id + 1) % n_verts], face[(v_id + 2) % n_verts] , face[(v_id + 3) % n_verts] });
                // Check for zero length edges and store them in a set
                const int current_edge_v1 = face[v_id];
                const int current_edge_v2 = face[(v_id + 1) % n_verts];
                if (fixed_vertices_set.find(current_edge_v1) != fixed_vertices_set.end() && fixed_vertices_set.find(current_edge_v2) != fixed_vertices_set.end()) continue; // Skip pinned edges
                TINYAD_ASSERT_NEQ(current_edge_v1, current_edge_v2);
                const double current_edge_length = (vertices[current_edge_v1] - vertices[current_edge_v2]).norm();
                mesh_edge_set.emplace(current_edge_v1, current_edge_v2, current_edge_length);
                if (current_edge_length <= 1e-16)
                {
                    TINYAD_INFO("Make Planar Faces detected an edge with near zero length! Merge vertices and try again.");
                    return vertices;
                }
                // Update min edge length
                min_edge_length = std::min(min_edge_length, current_edge_length);
            }
        }
        if (settings.verbose) TINYAD_INFO("Skipped " << skipped_faces << " fully pinned faces.")
        if (min_edge_length == std::numeric_limits<double>::infinity())
        {
            // This only happens when all vertices are pinned. We can return the input positions.
            return vertices;
        }
        // Compute a normalization factor to make the optimization invariant under mesh dimensions
        double normalization_factor = 1.0 / min_edge_length;
        if (settings.verbose) TINYAD_INFO("Make Planar Faces internal mesh scaling factor: " << normalization_factor);

        // Store all collected edges in a vector and apply the normalization factor + square lengths
        std::vector<MeshEdge> scaled_mesh_edges;
        for (const MeshEdge& e : mesh_edge_set)
        {
            scaled_mesh_edges.emplace_back(e.v1, e.v2, sqr(normalization_factor * e.length));
        }

        // Set up a TinyAD function
        int n = vertices.size();
        auto func = TinyAD::scalar_function<3>(TinyAD::range(n));

        // Compute base reduction matrix for fixed vertices
        // Matrix C maps from m-dim to n-dim space.
        // It is the identity map for all unconstrained vertices.
        SparseMatrix C;
        if (fixed_vertices.size() != 0)
        {
            const int m = 3 * (n - fixed_vertices.size());
            std::vector<bool> vertex_is_fixed_array = std::vector<bool>(n, false);
            for (int vertex_id : fixed_vertices)
            {
                vertex_is_fixed_array[vertex_id] = true;
            }

            std::vector<Triplet> C_triplets;
            C_triplets.reserve(m);
            int C_cols = 0;
            for (int vertex_id = 0; vertex_id < n; vertex_id++)
            {
                if (!vertex_is_fixed_array[vertex_id])
                {
                    C_triplets.emplace_back(3 * vertex_id + 0, C_cols++, 1.0);
                    C_triplets.emplace_back(3 * vertex_id + 1, C_cols++, 1.0);
                    C_triplets.emplace_back(3 * vertex_id + 2, C_cols++, 1.0);
                }
            }
            TINYAD_ASSERT_EQ(C_cols, m);

            C = SparseMatrix(3 * n, m);
            C.setFromTriplets(C_triplets.cbegin(), C_triplets.cend());
            TINYAD_ASSERT_EQ(C.rows(), 3 * n);
            TINYAD_ASSERT_EQ(C.cols(), m);
        }
        else
        {
            C = TinyAD::identity<double>(3 * n);
        }

        // Add closeness term
        double edge_objective_blend_factor = std::clamp(settings.edge_length_preservation_blend_factor, 0.0, 1.0);
        double closeness_weight = settings.initial_shape_preservation_weight;
        std::vector<Vec3d> normalized_vertices;
        normalized_vertices.reserve(n);
        for (const Vec3d& v : vertices)
        {
            normalized_vertices.push_back(normalization_factor * v);
        }
        std::vector<int> unconstrained_vertex_indices;
        for (int v_id = 0; v_id < n; v_id++)
        {
            if (fixed_vertices_set.find(v_id) == fixed_vertices_set.end())
            {
                unconstrained_vertex_indices.push_back(v_id);
            }
        }
        const int n_unconstrained_vertices = unconstrained_vertex_indices.size();
        if (settings.edge_length_preservation_blend_factor < 1 && n_unconstrained_vertices > 0)
        {
            func.add_elements<1>(TinyAD::range(n_unconstrained_vertices), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
            {
                // Evaluate element using either double or TinyAD::Double
                using T = TINYAD_SCALAR_TYPE(element);

                // Get the vertex position
                const int& current_v_index = unconstrained_vertex_indices[element.handle];
                Eigen::Vector3<T> current_v_pos = element.variables(current_v_index);
                Eigen::Vector3d orig_v_pos = normalized_vertices[current_v_index];

                // Compute squared distance
                return closeness_weight * (1.0 - edge_objective_blend_factor) * (current_v_pos - orig_v_pos).squaredNorm() / (double)n_unconstrained_vertices;
            });
        }

        // Add Edge length preservation + barrier term
        const int n_edges = scaled_mesh_edges.size();
        if (settings.edge_length_preservation_blend_factor > 0 && n_edges > 0)
        {
            func.add_elements<2>(TinyAD::range(n_edges), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
            {
                // Evaluate element using either double or TinyAD::Double
                using T = TINYAD_SCALAR_TYPE(element);

                // Get the edge
                const MeshEdge& current_edge = scaled_mesh_edges[element.handle];

                // Get the edge vertex positions
                Eigen::Vector3<T> v1_pos = element.variables(current_edge.v1);
                Eigen::Vector3<T> v2_pos = element.variables(current_edge.v2);
                
                T squared_edge_len = (v1_pos - v2_pos).squaredNorm();
                return closeness_weight * edge_objective_blend_factor * (squared_edge_len / current_edge.length + current_edge.length / squared_edge_len - 2.0) / (double)n_edges;
            });
        }

        // Add planarity term
        const int n_planarity_objectives = planarity_objective_indices.size();
        if (n_planarity_objectives > 0)
        {
            func.add_elements<4>(TinyAD::range(n_planarity_objectives), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
            {
                // Evaluate element using either double or TinyAD::Double
                using T = TINYAD_SCALAR_TYPE(element);

                const std::array<int, 4>& current_vertex_indices = planarity_objective_indices[element.handle];
                Eigen::Vector3<T> edge_a = (element.variables(current_vertex_indices[1]) - element.variables(current_vertex_indices[0])).normalized();
                Eigen::Vector3<T> edge_b = (element.variables(current_vertex_indices[2]) - element.variables(current_vertex_indices[1])).normalized();
                Eigen::Vector3<T> edge_c = (element.variables(current_vertex_indices[3]) - element.variables(current_vertex_indices[2])).normalized();

                return sqr(TinyAD::col_mat(edge_a, edge_b, edge_c).determinant()) / (double)n_planarity_objectives;
            });
        }

        // Init variables (scale the mesh to have a unit volume bounding box)
        Eigen::VectorXd x = func.x_from_data([&](long vertex_index)
        {
            return normalization_factor * vertices[vertex_index];
        });

        if (settings.verbose) TINYAD_INFO("Objective function setup done. Starting optimization...");

        // Optimize
        TinyAD::LinearSolver solver;
        double decay = 1.0;
        closeness_weight = settings.initial_shape_preservation_weight;
        if (settings.optimization_rounds > 1 && settings.initial_shape_preservation_weight != 0)
        {
            if (settings.target_shape_preservation_weight == 0)
            {
                // This rule is a bit random 
                decay = std::pow(std::min(0.1 * settings.initial_shape_preservation_weight, 1e-16) / settings.initial_shape_preservation_weight, 1.0 / (settings.optimization_rounds - 1.0));
            }
            else
            {
                decay = std::pow(settings.target_shape_preservation_weight / settings.initial_shape_preservation_weight, 1.0 / (settings.optimization_rounds - 1.0));
            }
        }
        for (int opt_round = 0; opt_round < settings.optimization_rounds; opt_round++)
        {
            if (settings.target_shape_preservation_weight == 0.0 && opt_round == settings.optimization_rounds - 1) closeness_weight = 0.0; 
            for (int iter = 0; iter < settings.max_iterations_per_round; iter++)
            {
                // eval function value, gradient and hessian
                auto [f, g, H_proj] = func.eval_with_hessian_proj(x, settings.projection_eps);

                // compute newton step direction
                Eigen::VectorXd d = TinyAD::newton_direction_reduced_basis(g, H_proj, C, solver, settings.w_identity);

                double newton_decrement = TinyAD::newton_decrement<double>(d, g);
                if (newton_decrement < settings.convergence_eps)
                {
                    if (settings.verbose) TINYAD_INFO("Newton decrement below convergence eps. Stopping early.");
                    break;
                }

                // line search for new x
                Eigen::VectorXd x_old = x;
                x = TinyAD::line_search(x, d, f, g, func, 1.0, 0.8, 128);
                
                // Early stopping 
                if (x == x_old)
                {
                    if (settings.verbose) TINYAD_INFO("Line search couldn't find improvement. Stopping early.");
                    break;
                }    
            }
            if (settings.verbose)
            {
                double f_val = func.eval(x);
                TINYAD_INFO("Energy | Shape Preservation Weight after round " << (opt_round + 1) << ": " << f_val << " | " << closeness_weight);
            } 
            // update closeness weight
            closeness_weight *= decay;
        }
        if (settings.verbose) TINYAD_INFO("Final energy: " << func.eval(x));

        // Extract solution (scale back to original dimensions)
        std::vector<Vec3d> optimized_vertex_positions = std::vector<Vec3d>(n, Vec3d::Zero());
        func.x_to_data(x, [&](int v_id, const Eigen::Vector3d& _p) {
            optimized_vertex_positions[v_id] = _p / normalization_factor;
        });

        return optimized_vertex_positions;
    }
}

PYBIND11_MODULE(mpfp, m)
{
    m.doc() = "Make Planar Faces Plus: A small geometry processing package for mesh planarization written in C++.";

    py::class_<MakePlanarFacesPlus::MakePlanarSettings>(m, "MakePlanarSettings")
        .def(py::init())
        .def_readwrite("optimization_rounds", &MakePlanarFacesPlus::MakePlanarSettings::optimization_rounds)
        .def_readwrite("max_iterations_per_round", &MakePlanarFacesPlus::MakePlanarSettings::max_iterations_per_round)
        .def_readwrite("initial_shape_preservation_weight", &MakePlanarFacesPlus::MakePlanarSettings::initial_shape_preservation_weight)
        .def_readwrite("target_shape_preservation_weight", &MakePlanarFacesPlus::MakePlanarSettings::target_shape_preservation_weight)
        .def_readwrite("edge_length_preservation_blend_factor", &MakePlanarFacesPlus::MakePlanarSettings::edge_length_preservation_blend_factor)
        .def_readwrite("verbose", &MakePlanarFacesPlus::MakePlanarSettings::verbose)
        .def_readwrite("projection_eps", &MakePlanarFacesPlus::MakePlanarSettings::projection_eps)
        .def_readwrite("w_identity", &MakePlanarFacesPlus::MakePlanarSettings::w_identity)
        .def_readwrite("convergence_eps", &MakePlanarFacesPlus::MakePlanarSettings::convergence_eps);
        
    m.def("make_planar_faces", &MakePlanarFacesPlus::make_planar_faces, "Continuous optimization that makes quad faces planar with minimal geometric loss.");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

}

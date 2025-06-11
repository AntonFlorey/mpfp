#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace MakePlanarFacesPlus
{

using Vec3d = Eigen::Vector3d;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

struct MeshEdge
{
	int v1, v2;
	double length;

	MeshEdge(int v1, int v2, double length)
	{
		this->v1 = v1 < v2 ? v1 : v2;
		this->v2 = v1 < v2 ? v2 : v1;
		this->length = length;
	}

	friend bool operator==(const MeshEdge& a, const MeshEdge& b) {
		return a.v1 == b.v1 && a.v2 == b.v2;
	}
};

inline void hash_combine(std::size_t& hash, const int& v);

struct MakePlanarSettings
{
	int optimization_rounds = 50;
	int max_iterations_per_round = 5;

	// Optimization settings
	double initial_shape_preservation_weight = 5.0;
	double target_shape_preservation_weight = 0.0;
	double edge_length_preservation_blend_factor = 0.5;

	// Optimizer settings
	bool verbose = true;
	double projection_eps = 1e-9;
	double w_identity = 1e-9;
	double convergence_eps = 1e-16;
};

std::vector<Vec3d> make_planar_faces(const std::vector<Vec3d>& vertices, const std::vector<std::vector<int>>& faces, const std::vector<int>& fixed_vertices, const MakePlanarSettings& settings);

}
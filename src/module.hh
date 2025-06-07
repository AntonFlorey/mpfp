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

struct MakePlanarSettings
{
	int optimization_rounds = 50;
	int max_iterations = 10;

	// Optimization settings
	double initial_closeness_weight = 5.0;
	double min_closeness_weight = 0.0;

	// Optimizer settings
	bool verbose = true;
	double projection_eps = 1e-9;
	double w_identity = 1e-9;
	double convergence_eps = 1e-16;
};

std::vector<Vec3d> make_planar_faces(const std::vector<Vec3d>& vertices, const std::vector<std::vector<int>>& faces, const std::vector<int>& fixed_vertices, const MakePlanarSettings& settings);

}
import mpfp
import pytest
import numpy as np

def test_version():
    assert hasattr(mpfp, "__version__")

def test_attributes():
    assert hasattr(mpfp, "make_planar_faces")
    assert hasattr(mpfp, "MakePlanarSettings")

def test_face_optimization():
    vertices = [[0,0,0], [1,0,0], [1,1,1], [0,1,0]]
    faces = [[0,1,2,3]]
    fixed_vertices = [0,1,3]
    opt_settings = mpfp.MakePlanarSettings()
    opt_settings.optimization_rounds = 100
    opt_settings.max_iterations_per_round = 5
    opt_settings.initial_shape_preservation_weight = 5
    opt_settings.target_shape_preservation_weight = 0
    opt_settings.edge_length_preservation_blend_factor = 0
    opt_settings.verbose = False
    opt_settings.projection_eps = 1e-9
    opt_settings.w_identity = 1e-9
    opt_settings.convergence_eps = 1e-16
    # optimize
    optimized_vertices = mpfp.make_planar_faces(vertices, faces, fixed_vertices, opt_settings)
    assert len(optimized_vertices) == 4
    assert np.allclose(optimized_vertices[2], [1,1,0])
    
if __name__ == "__main__":
    test_attributes()
    test_face_optimization()
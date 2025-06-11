# Make Planar Faces Plus
A small geometry processing package for mesh planarization written in C++.

<a href="https://github.com/patr-schm/TinyAD"><img src="https://img.shields.io/badge/Powered%20by-TinyAD-blue" alt="Badge referencing TinyAD." style="height:20px"/></a>

## Here is an example to get you started

```python
# necessary imports 
from mpfp import make_planar_faces, MakePlanarSettings
# prepare mesh data (in praxis you would create this data from your mesh)
vertices = [[0,0,0], [1,0,0], [1,1,0], [0,1,1]]
faces = [[0,1,2,3]]
fixed_vertices = [0,1,2]
# here is a list of all available settings (with default values):
opt_settings = MakePlanarSettings()
opt_settings.optimization_rounds = 50
opt_settings.max_iterations_per_round = 5
opt_settings.initial_shape_preservation_weight = 5
opt_settings.target_shape_preservation_weight = 0
opt_settings.edge_length_preservation_blend_factor = 0.5
opt_settings.verbose = True
opt_settings.projection_eps = 1e-9
opt_settings.w_identity = 1e-9
opt_settings.convergence_eps = 1e-16
# optimize
optimized_vertices = make_planar_faces(vertices, faces, fixed_vertices, opt_settings)
# print the result
print(optimized_vertices)
```

## How to encode your Mesh
You provide your mesh to the `make_planar_faces` function via the two parameters:

- `vertices`: A list of 3D vertex coordinates. You can provide them as a 2d list or a numpy array with shape (n, 3).
- `faces`: A list of all mesh faces. Each face has to be provided as a list of vertex indices in ccw or cw order.

## Details
The function provided by this module aims to make each face of a mesh planar. It solves a global optimization problem in order to make faces planar while preserving the objects shape as much as possible.

You can control the strength of this shape preservation objective via the `initial_shape_preservation_weight` and `target_shape_preservation_weight` parameter. The algorithm will interpolate between the two while optimizing. 

There are two objective terms that both favor shape preservation. One penalizes the distance of all vertices to their original position and the other one penalizes changes in edge lengths. You can blend these two objectives via the `edge_length_preservation_blend_factor` parameter: A value of 0 will turn off the edge length objective, a value of 1 will turn off the vertex distance objective. The default value 0.5 will lead to a balance between the two.

How fast the shape preservation weight decays is controlled by the `optimization_rounds` setting. If set to 2, the first round will use the initial weight, the second round the target weight. More rounds will add more steps inbetween these two values, leading to a more graceful descent. 

The `max_iterations_per_round` setting determines how many optimization steps are performed per round (so for one weight value). A round will be stopped early if the objective funtion improvement falls below the `convergence_eps` threshold. 

Here is an overview of the optimization process:

```python
# python pseudo-code of the optimization process
shape_weight = "initial_shape_preservation_weight"
for opt_round in range("optimization_rounds"):
    for opt_step in range("max_iterations_per_round"):
        do_optimization_step()
        if improvement < "convergence_eps":
            break
    # the decay factor is chosen such that
    # shape_weight = "target_shape_preservation_weight" in the last round
    shape_weight *= decay_factor
```

The algorithm will always try to optimize the entire mesh. By providing a list of pinned vertx indices, all selected vertices will not be affected by the operator. This may be useful when you want to preserve certain features.
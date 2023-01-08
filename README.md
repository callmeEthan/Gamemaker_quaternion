# Quaternion for gamemaker
These are some useful [quaternion function](https://github.com/callmeEthan/Gamemaker_quaternion/blob/main/scripts/Quaternion/Quaternion.gml) for gamemaker to use in place of GM traditional Euler angles.  
Quaternion can combine multiple rotation without encountering Gimbal Lock.  
The quaternion is format as [x, y, z, w], this is how you usually supply to GPU for vertex transform.  

# Transform matrix
You can create a transform matrix from a quaternion that you can use for rendering.
```
matrix = matrix_build_quaternion(x, y, z, quaternion, xscale, yscale, zscale);
matrix_set(matrix_world, matrix);
vertex_submit(model, pr_trianglelist, texture);
```
# Vertex transform
There is also a [shader](https://github.com/callmeEthan/Gamemaker_quaternion/blob/main/shaders/sh_quat_transform/sh_quat_transform.vsh) for vertex transform, if you'd like that route instead.

```
shader_set(sh_quat_transform)
shader_set_uniform_f(shader_get_uniform(sh_quat_transform, "u_position"), x, y, z);
shader_set_uniform_f(shader_get_uniform(sh_quat_transform, "u_scale"), xscale, yscale, zscale);
shader_set_uniform_f(shader_get_uniform(sh_quat_transform, "u_rotate"), q[0], q[1], q[2], q[3]);
vertex_submit(model, pr_trianglelist, texture);
```

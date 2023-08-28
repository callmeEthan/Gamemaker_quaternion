# Quaternion for gamemaker
These are some useful [quaternion function](https://github.com/callmeEthan/Gamemaker_quaternion/blob/main/scripts/Quaternion/Quaternion.gml) for gamemaker to use in place of GM traditional Euler angles.  
Quaternion can combine multiple rotation without suffering from Gimbal Lock, which is a well-known problem for euler angles.  

# Replacing Euler angles
**Basic quaternion**  
The quaternion is format as [x, y, z, w] (1d array), this is how you usually supply to GPU for vertex transform. Default quaternion with no rotation is [0, 0, 0, 1].   
```
rotation = quaternion_identity();
```
**Quaternion rotation**  
You can rotate a quaternion by multiplying two quaternions.
```
// Rotate a quaternion around world Z-axis
var quat = angle_to_quaternion(0, 0, angle);
rotation = quaternion_multiply(quat, rotation);
```
You can also rotate a quaternion by it's own local axis.  
```
// Rotate a quaternion around local Y-axis
var quat = angle_to_quaternion(0, angle, 0);
rotation = quaternion_multiply(rotation, quat);
```


# Transform matrix
You can create a matrix from a quaternion that you can use for transform.
```
matrix = matrix_build_quaternion(x, y, z, quaternion, xscale, yscale, zscale);
matrix_set(matrix_world, matrix);
vertex_submit(model, pr_trianglelist, texture);
```
# Vertex transform
There is also a [shader](https://github.com/callmeEthan/Gamemaker_quaternion/blob/main/shaders/sh_quat_transform/sh_quat_transform.vsh) for vertex transform, if you like that route instead.

```
shader_set(sh_quat_transform);
shader_set_uniform_f(shader_get_uniform(sh_quat_transform, "u_position"), x, y, z);
shader_set_uniform_f(shader_get_uniform(sh_quat_transform, "u_scale"), xscale, yscale, zscale);
shader_set_uniform_f_array(shader_get_uniform(sh_quat_transform, "u_rotate"), quaternion);
vertex_submit(model, pr_trianglelist, texture);
shader_reset();
```

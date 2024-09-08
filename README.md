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
**Interpolation**  
You can interpolate between two quaternion unit using one of the three function
```
var q1, q2; 				// Two quaternion rotation unit
var amount = 2 				// interpolation value between 0 and 1
var q = array_create(4); 	// Final output
quaternion_lerp(q1, q2, amount, q);		// linear interpolation
quaternion_slerp(q1, q2, amount, q);	// spherical-linear-interpolation
quaternion_nlerp(q1, q2, amount, q);	// normalized-linear-interpolation.
```
- Linear interpolation is fast, but does not represent angular rotation correctly (low quality).  
- Spherical linear interpolation produce correct angular rotation with constant velocity (high quality).  
- Normalized linear interpolation is faster than slerp, correctly represent angular rotation but does not perform a constant velocity interpolation (variable acceleration/decceleration during interpolation).  

**Transformation**  
Transform a unit vector by a quaternion
```
var pos = quaternion_transform_vector(quat, x, y, z)
```
Result is an array containing new coordinate of the unit vector [x,y,z].  
**Difference**  
You can get the rotation difference between two quaternion unit:
```
var q = quaternion_difference(q1, q2);	// Return array quaternion unit [x,y,z,w]
```
Or get quaternion unit between two direction vector:
```
var q = quaternion_vector_angle(v1, v2);   // Return array quaternion unit [x,y,z,w]
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

/// Quaternion scripts
/// 08/09/2023
/// @callmeEthan

function quaternion_identity(array=array_create(4))
{
	array[@0]=0
	array[@1]=0
	array[@2]=0
	array[@3]=1.0001
	return array;
}

function quaternion_build(xaxis, yaxis, zaxis, angle, array=array_create(4))
{
	gml_pragma("forceinline");
	angle = degtorad(angle)
	var s = sin(angle/2);
	array[@0] = xaxis * s;
	array[@1] = yaxis * s;
	array[@2] = zaxis * s;
	array[@3] = cos(angle/2);
	return array;
}

function quaternion_to_angle(q) {
	gml_pragma("forceinline");
    // roll (x-axis rotation)
    var sinr_cosp = 2 * (q[3] * q[0] + q[1] * q[2]);
    var cosr_cosp = 1 - 2 * (q[0] * q[0] + q[1] * q[1]);
    var roll = arctan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    var sinp = sqrt(1 + 2 * (q[3] * q[0] - q[1] * q[2]));
    var cosp = sqrt(1 - 2 * (q[3] * q[0] - q[1] * q[2]));
    var pitch = 2 * arctan2(sinp, cosp) - pi / 2;

    // yaw (z-axis rotation)
    var siny_cosp = 2 * (q[3] * q[2] + q[0] * q[1]);
    var cosy_cosp = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
    var yaw = arctan2(siny_cosp, cosy_cosp);

	yaw = radtodeg(yaw);
	pitch = radtodeg(pitch);
	roll = radtodeg(roll);
    return [roll, pitch, yaw];
}

function angle_to_quaternion(xangle, yangle, zangle, output=array_create(4))
{
	gml_pragma("forceinline");
    // Abbreviations for the various angular functions
    var cr = dcos(xangle * 0.5);
    var sr = dsin(xangle * 0.5);
    var cp = dcos(yangle * 0.5);
    var sp = dsin(yangle * 0.5);
    var cy = dcos(zangle * 0.5);
    var sy = dsin(zangle * 0.5);
	
	var q = output;
    q[@3] = cr * cp * cy + sr * sp * sy;
    q[@0] = sr * cp * cy - cr * sp * sy;
    q[@1] = cr * sp * cy + sr * cp * sy;
    q[@2] = cr * cp * sy - sr * sp * cy;
    return q;
}

function axis_to_quaternion(xaxis, yaxis, zaxis, angle)
{
	return quaternion_build(xaxis, yaxis, zaxis, angle)
}

function quaternion_multiply(R, S, array = array_create(4))
{ 
	gml_pragma("forceinline");
	var Qx = R[3] * S[0] + R[0] * S[3] + R[1] * S[2] - R[2] * S[1];
	var Qy = R[3] * S[1] + R[1] * S[3] + R[2] * S[0] - R[0] * S[2];
	var Qz = R[3] * S[2] + R[2] * S[3] + R[0] * S[1] - R[1] * S[0];
	var Qw = R[3] * S[3] - R[0] * S[0] - R[1] * S[1] - R[2] * S[2];
	array[@0] = Qx;
	array[@1] = Qy;
	array[@2] = Qz;
	array[@3] = Qw;
	return array;
}

function quaternion_rotate_local(q, xrot, yrot, zrot)
{
	// Rotate a quaternion around it's local axis.
	var rot = angle_to_quaternion(xrot, yrot, zrot)
	return quaternion_multiply(q, rot, q);
}

function quaternion_rotate_world(q, xrot, yrot, zrot)
{
	// Rotate a quaternion around world axis.
	var rot = angle_to_quaternion(xrot, yrot, zrot)
	return quaternion_multiply(rot, q, q);
}

function matrix_build_quaternion(x, y, z, quaternion, xscale, yscale, zscale, matrix=array_create(16))
{
	// Build transform matrix based on quaternion rotation instead of Euler angle
	// https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
	//	You can defined output array (1D array length 16), it will write transform matrix directly onto output array and return nothing.
	//	Output should be an array size 16 to write data directly, if not defined return a new array.
   gml_pragma("forceinline");
   var mat = matrix;
   var sqw = quaternion[3]*quaternion[3];
   var sqx = quaternion[0]*quaternion[0];
   var sqy = quaternion[1]*quaternion[1];
   var sqz = quaternion[2]*quaternion[2];
   mat[@ 0] = (sqx - sqy - sqz + sqw) * xscale; // since sqw + sqx + sqy + sqz =1
   mat[@ 5] = (-sqx + sqy - sqz + sqw) * yscale;
   mat[@ 10] = (-sqx - sqy + sqz + sqw) * zscale;
   
   var tmp1 = quaternion[0]*quaternion[1];
   var tmp2 = quaternion[2]*quaternion[3];
   mat[@ 1] = 2.0 * (tmp1 + tmp2) * xscale;
   mat[@ 4] = 2.0 * (tmp1 - tmp2) * yscale;
   
   tmp1 = quaternion[0]*quaternion[2];
   tmp2 = quaternion[1]*quaternion[3];
   mat[@ 2] = 2.0 * (tmp1 - tmp2) * xscale;
   mat[@ 8] = 2.0 * (tmp1 + tmp2) * zscale;
   
   tmp1 = quaternion[1]*quaternion[2];
   tmp2 = quaternion[0]*quaternion[3];
   mat[@ 6] = 2.0 * (tmp1 + tmp2) * yscale;
   mat[@ 9] = 2.0 * (tmp1 - tmp2) * zscale;
	
	mat[@ 12] = x;
	mat[@ 13] = y;
	mat[@ 14] = z;
	mat[@ 15] = 1.0;
	return mat;
}

function quaternion_conjugate(q) 
{
	//	Return the conjugate of a quaternion 
	q[@0] *= -1;
	q[@1] *= -1;
	q[@2] *= -1;
	return q;
}

function quaternion_to_rotation_matrix(q)
{
	// Convert quaternion rotate to rotation matrix
	//https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
	var M = array_create(9)
    var sqw = q[3]*q[3];
    var sqx = q[0]*q[0];
    var sqy = q[1]*q[1];
    var sqz = q[2]*q[2];

    // invs (inverse square length) is only required if quaternion is not already normalised
    var invs = 1 / (sqx + sqy + sqz + sqw)
    M[@0] = ( sqx - sqy - sqz + sqw)*invs ; // since sqw + sqx + sqy + sqz =1/invs*invs
    M[@4] = (-sqx + sqy - sqz + sqw)*invs ;
    M[@8] = (-sqx - sqy + sqz + sqw)*invs ;
    
    var tmp1 = q[0]*q[1];
    var tmp2 = q[2]*q[3];
    M[@3] = 2.0 * (tmp1 + tmp2)*invs ;
    M[@1] = 2.0 * (tmp1 - tmp2)*invs ;
    
    tmp1 = q[0]*q[2];
    tmp2 = q[1]*q[3];
    M[@6] = 2.0 * (tmp1 - tmp2)*invs ;
    M[@2] = 2.0 * (tmp1 + tmp2)*invs ;
    tmp1 = q[1]*q[2];
    tmp2 = q[0]*q[3];
    M[@7] = 2.0 * (tmp1 + tmp2)*invs ;
    M[@5] = 2.0 * (tmp1 - tmp2)*invs ;      
	
	return M
}

function quaternion_normalize(q, array=array_create(4))
{
	// Normalize vec4 quaternion, output will be written directly onto result array[4].
	//	output result can be the same as input, if result array not defined, return a new array.
	gml_pragma("forceinline");
	var l = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
	if l == 0 {return quaternion_identity()}
	l = sqrt(l);
	var j = 1 / l;
	array[@0] = q[0] * j;
	array[@1] = q[1] * j;
	array[@2] = q[2] * j;
	array[@3] = q[3] * j;
	return array;
}

function quaternion_lerp(qa, qb, amount, array=array_create(4))
{
	// Linear interpolation between two vector, fast but lower quality (use slerp for better quality).
	array[@0] = lerp(qa, qb, amount);
	array[@1] = lerp(qa, qb, amount);
	array[@2] = lerp(qa, qb, amount);
	array[@3] = lerp(qa, qb, amount);;
	return array;
}

function quaternion_slerp(qa, qb, amount, array=array_create(4))
{
	// Interpolate between two quaternions representing rotations with Spherical-linear-interpolation for better quality
	var qax=qa[0], qay=qa[1], qaz=qa[2], qaw=qa[3];
	var qbx=qb[0], qby=qb[1], qbz=qb[2], qbw=qb[3];
	var dot = qaw * qbw + qax * qbx + qay * qby + qaz * qaz;
	if dot>0.99
	{	// Linear interpolation
		array[@0] = lerp(qax, qbx, amount);
		array[@1] = lerp(qay, qby, amount);
		array[@2] = lerp(qaz, qbz, amount);
		array[@3] = lerp(qaw, qbw, amount);
		return array;
	}
	var angle = arccos(dot);
	var denom = sin(angle);
	// Spherical linear interpolation.
	var r1 = sin((1-amount)*angle);
	var r2 = sin(amount*angle);
	array[@0] = (qax * r1 + qbx * r2)/denom
	array[@1] = (qay * r1 + qby * r2)/denom
	array[@2] = (qaz * r1 + qbz * r2)/denom
	array[@3] = (qaw * r1 + qbw * r2)/denom
	return array;
}

function quaternion_nlerp(qa, qb, amount, array)
{
	// Interpolate between two quaternions representing rotations with Normalized-linear-interpolation.
	// While nlerp is faster and more efficient than slerp; however, it does not perform a constant angular velocity interpolation (there might be variable acceleration/decceleration during interpolation)
	// It should be fine for short interpolation intervals (small angles).
	var qx = lerp(qa[0], qb[0], amount);
	var qy = lerp(qa[1], qb[1], amount);
	var qz = lerp(qa[2], qb[2], amount);
	var qw = lerp(qa[3], qb[3], amount);
	var l = sqrt(qx*qx + qy*qy + qz*qz + qw*qw);
	array[@0] = qx / l;
	array[@1] = qy / l;
	array[@2] = qz / l;
	array[@3] = qw / l;
	return array;
}

function quaternion_transform_vector(quat, x, y, z, array=array_create(3))
{
	// Transform a vector by quaternion unit
	var qx = quat[0]
	var qy = quat[1]
	var qz = quat[2]
	var qw = quat[3]
	array[@0] = qw*qw*x + 2*qy*qw*z - 2*qz*qw*y + qx*qx*x + 2*qy*qx*y + 2*qz*qx*z - qz*qz*x - qy*qy*x;
	array[@1] = 2*qx*qy*x + qy*qy*y + 2*qz*qy*z + 2*qw*qz*x - qz*qz*y + qw*qw*y - 2*qx*qw*z - qx*qx*y;
	array[@2] = 2*qx*qz*x + 2*qy*qz*y + qz*qz*z - 2*qw*qy*x - qy*qy*z + 2*qw*qx*y - qx*qx*z + qw*qw*z;
	return array;
}

function quaternion_vector_angle(v0, v1, array=array_create(4))
{
	// Return rotation unit between two directional vector (nearest rotation angle).
	// v1,v2: is array[3] represent 3d direction of the vectors [x,y,z].
	// array: result array to write quaternion angle; if no array is defined return new array.
	var d = dot_product_3d(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2]);
	var m1 = sqrt(power(v0[0],2) + power(v0[1],2) + power(v0[2],2));
	var m2 = sqrt(power(v1[0],2) + power(v1[1],2) + power(v1[2],2));
	var angle = arccos( d/(m1*m2)); 
	
	if angle==0 {quaternion_identity(array); return array};
	if angle==pi/2 {if abs(v0[2])=1 cross_product(v0, [0,1,0], array) else cross_product(v0, [0,0,1], array)
	} else {cross_product(v0, v1, array);}
	normalize(array);
	
	var s = sin(angle/2);
	array[@0]*= s;
	array[@1]*= s;
	array[@2]*= s;
	array[@3] = cos(angle/2);
	return array;
}

function quaternion_difference(qa, qb, array=array_create(4))
{
	// Difference of 2 quaternion unit, like angle difference
	quaternion_conjugate(qb);
	quaternion_multiply(qb, qa, array);
	if qb!=array quaternion_conjugate(qb);
}

//
// Simple passthrough vertex shader
//
attribute vec3 in_Position;                  // (x,y,z)
attribute vec3 in_Normal;                  // (x,y,z)     unused in this shader.
attribute vec4 in_Colour;                    // (r,g,b,a)
attribute vec2 in_TextureCoord;              // (u,v)

varying vec2 v_vTexcoord;
varying vec4 v_vColour;
varying vec3 v_vNormal;

uniform vec3 u_position;
uniform vec3 u_scale;
uniform vec4 u_rotate;

vec3 rotate_quaternion(vec3 vec, vec4 q)	// ROTATE BY QUATERNION 
{
	return vec + 2.0*cross(cross(vec, q.xyz ) + q.w*vec, q.xyz);
}
vec3 transform_vertex(vec3 vec, vec3 pos, vec4 rot, vec3 scale)
{
	vec3 vertex = vec * scale;
	vertex = rotate_quaternion(vertex, rot);
	return vertex + pos;
}

void main()
{
	
	vec3 object_space_pos = transform_vertex(in_Position, u_position, u_rotate, u_scale);
    gl_Position = gm_Matrices[MATRIX_WORLD_VIEW_PROJECTION] * vec4(object_space_pos, 1.0);
    
    v_vColour = in_Colour;
    v_vTexcoord = in_TextureCoord;
	v_vNormal = in_Normal;
}

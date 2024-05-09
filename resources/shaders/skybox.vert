#version 330 core

layout(location = 0) in vec3 pos; // Position of the vertex
//layout(location = 1) in vec3 normal;   // Normal of the vertex
//layout(location = 3) in vec3 texCoords;

uniform mat4 view, projection, rotation;

out vec3 tex_coord;



void main() {

    tex_coord = vec3(pos.x, pos.y, -pos.z);
    vec4 world_pos = projection*view*rotation*vec4(pos, 1.0);
    gl_Position = vec4(world_pos.x, world_pos.y, world_pos.w, world_pos.w);

//    tex_coord = pos;
  // gl_Position = projection*view*vec4(pos, 1.0);
}

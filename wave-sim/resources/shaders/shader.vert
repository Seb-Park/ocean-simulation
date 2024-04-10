#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex
layout(location = 2) in vec3 force;

uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;

uniform mat3 inverseTransposeModel;

out vec4 normal_worldSpace;
out vec4 position_worldSpace;
out vec4 force_worldSpace;

void main() {
    normal_worldSpace   = vec4(normalize(inverseTransposeModel * normal), 0);
    position_worldSpace = vec4(position, 1.0);
    force_worldSpace = vec4(force, 1.0);
//    force_worldSpace = vec4(1.0, 0.0, 0.0, 1.0);

    gl_Position = proj * view * model * vec4(position, 1.0);
}

#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex

uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;
uniform mat4 inverseView;

uniform mat3 inverseTransposeModel;

out vec3 normal_cameraSpace;
out vec3 normal_worldSpace;
out vec3 camera_worldSpace;
out vec3 pos;

void main() {
    normal_cameraSpace = normalize(inverse(transpose(mat3(view))) * inverseTransposeModel * normal);
    camera_worldSpace = vec3(inverseView * vec4(0.f, 0.f, 0.f, 1.f));
    normal_worldSpace = normal;
    pos = position;

    gl_Position = proj * view * model * vec4(position, 1);
}

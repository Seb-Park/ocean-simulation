#version 330 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
out vec3 col;

void main() {
   gl_Position = vec4(position, 1.0);
   col = normal;
}

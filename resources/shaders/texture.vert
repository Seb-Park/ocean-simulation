#version 330 core

// TASK 15: add a second layout variable representing a UV coordinate
layout (location = 0) in vec3 position;
layout (location = 1) in vec2 uv_raw;

// TASK 16: create an "out" variable representing a UV coordinate
out vec2 uv;

void main() {
    // TASK 16: assign the UV layout variable to the UV "out" variable
    uv = uv_raw;

    gl_Position = vec4(position, 1.0);
}

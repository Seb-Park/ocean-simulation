#version 330 core

in vec3 oldPos;
in vec3 newPos;
in vec3 ray;
out vec4 fragColor;

void main() {
    float oldArea = length(dFdx(oldPos)) * length(dFdy(oldPos));
    float newArea = length(dFdx(newPos)) * length(dFdy(newPos));
    fragColor = vec4(oldArea / newArea * 0.2, 1.0, 0.0, 0.0);
}

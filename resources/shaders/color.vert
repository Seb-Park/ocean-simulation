#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex
layout(location = 3) in vec3 texCoords;   // Normal of the vertex

uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;
uniform mat4 inverseView;
uniform vec3 light = vec3(0.0, -1.0, 0.0);

const float IOR_AIR = 1.0;
const float IOR_WATER = 1.333;

uniform sampler2D water;

out vec3 col;
out vec3 oldPos;
out vec3 newPos;
out vec3 ray;

vec3 project(vec3 origin, vec3 ray, vec3 refractedLight) {
    origin += ray;
    float tplane = (-origin.y - 1.0) / refractedLight.y;
    return origin + refractedLight * tplane;
}

void main() {
    vec4 info = texture(water, position.xy * 0.5 + 0.5);
    info.ba *= 0.5;
    vec3 normal = vec3(info.b, sqrt(1.0 - dot(info.ba, info.ba)), info.a);

    /* project the vertices along the refracted vertex ray */
    vec3 refractedLight = refract(-light, vec3(0.0, 1.0, 0.0), IOR_AIR / IOR_WATER);
    ray = refract(-light, normal, IOR_AIR / IOR_WATER);
    oldPos = project(position.xzy, refractedLight, refractedLight);
    newPos = project(position.xzy + vec3(0.0, info.r, 0.0), ray, refractedLight);

    gl_Position = vec4(0.75 * (newPos.xz + refractedLight.xz / refractedLight.y), 0.0, 1.0);
}

#version 330 core

// All in object space
layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 norm;
layout (location = 2) in vec2 uv;
layout (location = 3) in vec3 color;

uniform mat4 model, view, projection;
uniform vec3 worldSpace_camPos;

//out vec3 worldSpace_norm;
out vec2 tex_coord;
out vec3 vertColor;
out float visibility;
out vec4 clipSpace;
out vec3 toCameraVector;

// fog
const float density = .01f;
const float gradient = 4.f;

const float tiling = 6.f;

vec4 worldSpace_pos;

// water referenced from: ThinMatrix https://www.youtube.com/watch?v=qgDPSnZPGMA

void main() {
    worldSpace_pos = view*model*vec4(pos, 1.0);
    //worldSpace_norm = vec3(transpose(inverse(model))*vec4(norm, 0.0));

    vec4 positionRelationToCam = view * vec4(pos, 1.f);
    float distance = length(positionRelationToCam.xyz);
    visibility = exp(-pow((distance*density), gradient));
    visibility = clamp(visibility, 0.f, 1.f);

    tex_coord = uv*tiling;
    vertColor = color;

    clipSpace = projection*worldSpace_pos;
    gl_Position = clipSpace;

    toCameraVector = worldSpace_camPos - worldSpace_pos.xyz;
}

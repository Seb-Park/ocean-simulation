#version 330 core

// All in object space
layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 norm;
layout (location = 2) in vec2 uv;
layout (location = 3) in vec3 color;

uniform mat4 model, view, projection;
uniform vec3 worldSpace_camPos;

out vec3 worldSpace_pos;
out vec3 worldSpace_norm;
out vec2 tex_coord;
out vec3 vertColor;
out float visibility;

// fog
const float density = .01f;
const float gradient = 4.f;

// REFLECTION/REFRACTION
uniform vec4 plane; // horizontal plane

void main() {
    worldSpace_pos = vec3(view*model*vec4(pos, 1.0));
    worldSpace_norm = vec3(transpose(inverse(model))*vec4(norm, 0.0));

    vec4 positionRelationToCam = view * vec4(pos, 1.f);
    float distance = length(positionRelationToCam.xyz);
    visibility = exp(-pow((distance*density), gradient));
    visibility = clamp(visibility, 0.f, 1.f);

   // gl_ClipDistance[0]=dot(vec4(worldSpace_pos, 1.f), plane);


    tex_coord = uv;
    vertColor = color;

    gl_Position = projection*view*model*vec4(pos, 1.0);
}

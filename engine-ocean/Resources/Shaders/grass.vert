#version 330 core
layout (location = 0) in vec3 pos;

out VS_OUT {
  float visibility;
} vs_out;

uniform mat4 view, projection;
uniform vec4 plane;

// fog
const float density = .01f;
const float gradient = 4.f;

void main() {

    gl_ClipDistance[0]=dot(vec4(pos, 1.f), plane);

    vec4 positionRelationToCam = view * vec4(pos, 1.f);
    float distance = length(positionRelationToCam.xyz);
    float vis = exp(-pow((distance*density), gradient));
    vs_out.visibility = clamp(vis, 0.f, 1.f);

    gl_Position = vec4(pos, 1.0);
}

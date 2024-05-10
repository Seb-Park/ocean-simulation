#version 410 core
out vec4 fragColor;

in vec3 normal_worldSpace;
in vec3 pos;
in vec3 newPos;
in vec4 col;
in float refractProb;

uniform sampler2D normSamp;

void main() {
//    fragColor = vec4(vec3((pos[0] + 1) / 2, (pos[1] + 1) / 2, 0.f), 1.f);
    float oldArea = length(dFdx(vec3(pos[0], pos[2], pos[1]))) * length(dFdy(vec3(pos[0], pos[2], pos[1])));
    float newArea = length(dFdx(vec3(newPos[0], newPos[2], newPos[1]))) * length(dFdy(vec3(newPos[0], newPos[2], newPos[1])));
    float areaRatio = oldArea / newArea;
    float intensity = pow(areaRatio * .3f, 1.5f);
    fragColor = vec4(0.98, 1, .78, intensity * refractProb);
//    fragColor = col;
}

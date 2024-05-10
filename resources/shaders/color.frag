#version 410 core
out vec4 fragColor;

in vec3 normal_cameraSpace;
in vec3 camera_worldSpace;
in vec3 normal_worldSpace;
in vec3 pos;
in vec3 oldPosFlat;
in float intensity;

//uniform float multiplier = 1.f;
//uniform float contrast = 10.f;
//uniform float intExp = 0.5f;
//uniform float scale = 1.f;

//uniform float multiplier = .9f;
//uniform float contrast = 20.f;
//uniform float intExp = 0.f;
//uniform float scale = 1.f;

uniform float multiplier = .5f;
uniform float contrast = 1.5f;
uniform float intExp = 0.f;
uniform float scale = 1.f;

//uniform vec4 baseColor =vec4();

void main() {
    fragColor = 0.5f * (vec4(normal_worldSpace[0], 0.f, normal_worldSpace[2], 1.f) + 1);
    float oldArea = length(dFdx(oldPosFlat)) * length(dFdy(oldPosFlat));
    float newArea = length(dFdx(pos)) * length(dFdy(pos));
    float areaRatio = oldArea / newArea * multiplier;
//    if(oldArea / newArea > 3) {
//        areaRatio = 0.f;
//    }
//    fragColor = clamp(vec4(intensity, intensity, (1.f - intensity) * .5f, intensity), 0.f, 1.f);
    fragColor = vec4(pow(areaRatio, contrast)) * pow(intensity, intExp);
    float finalInt = pow(areaRatio, contrast) * pow(intensity, intExp);
    finalInt = clamp(finalInt, 0, 1);
    fragColor = vec4(1, 1, 1, finalInt * scale);
//    fragColor = vec4(vec3(1), pow(oldArea / newArea * .2f, 1.5f));
    fragColor = vec4((normal_worldSpace + 1) / 2, 1.f);
}

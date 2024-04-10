#version 330 core
out vec4 fragColor;

in vec3 normal_cameraSpace;
in vec3 camera_worldSpace;
in vec3 normal_worldSpace;
in vec3 pos;

uniform int   wire  = 0;
uniform float red   = 1.0;
uniform float green = 1.0;
uniform float blue  = 1.0;
uniform float alpha = 1.0;

void main() {
    // Do lighting in camera space
    vec3 lightDir = normalize(vec3(0, 0.5, 1));
    float d = clamp(dot(normal_cameraSpace, lightDir), 0, 1);
    vec3 reflectedLight = lightDir - 2 * dot(lightDir, normal_worldSpace) * normal_worldSpace;
    vec3 posToCam = normalize(camera_worldSpace - pos);
    float spec = pow(dot(posToCam, reflectedLight), 2.f);

    fragColor = clamp(0.5f * vec4(red * d, green * d, blue * d, 0.5f) + 0.5f * vec4(1, 1, 1, 1) * spec, 0, 1);
}

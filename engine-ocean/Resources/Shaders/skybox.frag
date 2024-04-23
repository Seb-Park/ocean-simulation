#version 330 core

in vec3 tex_coord;

uniform samplerCube cubeMap;
uniform vec3 skyColor;

out vec4 fragColor;

const float lowerLimit = -50.f;
const float upperLimit = 100.f;


void main() {

    //fragColor = vec4(1.f);
    vec4 finalColor = texture(cubeMap, tex_coord);

    // blending bottom of skybox to skyColor
    float factor = (tex_coord.y - lowerLimit) / (upperLimit - lowerLimit);
    factor = clamp(factor, 0.f, 1.f);
    fragColor = mix(vec4(skyColor, 1.f), finalColor, factor);
}

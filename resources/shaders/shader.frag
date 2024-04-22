#version 330 core
out vec4 fragColor;

in vec3 normal_cameraSpace;
in vec3 camera_worldSpace;
in vec3 normal_worldSpace;
in vec3 pos;
in vec3 refrPos;
in float refrProb;
in vec2 uv;

uniform int   wire  = 0;
uniform float red   = 1.0;
uniform float green = 1.0;
uniform float blue  = 1.0;
uniform float alpha = 1.0;
uniform sampler2D sampler;
uniform vec2 widthBounds;
uniform vec2 lengthBounds;
//uniform float test = 0;

vec2 uvFromWorldPoint(vec3 point) {
    float u = (point.x - widthBounds[0]) / (widthBounds[1] - widthBounds[0]);
    float v = (point.z - lengthBounds[0]) / (lengthBounds[1] - lengthBounds[0]);
    return vec2(u, v);
}

void main() {
    // Do lighting in camera space
    vec3 lightDir = normalize(vec3(0, 0.5, 1));
    lightDir = normalize(vec3(0.f, 3.f, 0.f) - pos);
//    float d = clamp(dot(normal_cameraSpace, lightDir), 0, 1);
    float d = clamp(dot(normal_worldSpace, lightDir), 0, 1);
    vec3 reflectedLight = lightDir - 2 * dot(lightDir, normal_worldSpace) * normal_worldSpace;
    vec3 posToCam = normalize(camera_worldSpace - pos);
    float spec = pow(dot(posToCam, reflectedLight), 2.f);

//    fragColor = texture(sampler, vec2(0.5f, 0.5f));
//    fragColor = vec4(abs(pos.x / 160.f), pos.y, 0.f, 1.f);
//    fragColor = vec4(uv.y, uv.y, 0.f, 1.f);
//    fragColor = vec4(camera_worldSpace.x - pos[0], camera_worldSpace.y - pos[1], pos[2], 1.f);
//    fragColor = vec4(- pos[0], 0.f, 0.f, 1.f);
//    fragColor = vec4((pos - vec3(widthBounds[0], 0, lengthBounds[0])) / 5.f, 1.f);
//    fragColor = vec4(fragColor.x, 0.f, fragColor.z, 1.f);
//    fragColor = vec4(test, test, test, 1.f);
    vec2 refrUV = uvFromWorldPoint(refrPos);
    vec4 transmissive = vec4(vec3(refrUV, 1.f - refrUV.y), 1.f);

    fragColor = 0.25f * vec4(red * d, green * d, blue * d, 1.0f); // Diffuse
    fragColor += 0.75f * vec4(1, 1, 1, 1) * pow(spec, 10.f); // Specular TODO: Pass multiplications as uniforms.
    fragColor = clamp(fragColor, 0.f, 1.f); // Clamp
    fragColor *=  (1 - (refrProb / 1.f));
    fragColor += (refrProb / 1.5f) * transmissive;
//    fragColor = transmissive * refrProb;
    fragColor = vec4(vec3(fragColor), 1.5f);
    // Dividing refrProb by 2 just for heuristic. Want more phong to show through.
//    fragColor = clamp(fragColor, 0.f, 1.f);
//    fragColor = vec4(refrProb, 0.f, 0.f, 1.f);
}

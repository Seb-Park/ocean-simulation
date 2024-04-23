#version 330 core
out vec4 fragColor;

in vec3 normal_cameraSpace;
in vec3 camera_worldSpace;
in vec3 normal_worldSpace;
in vec3 pos;
in vec3 refrPos;
in float refrProb;
in vec2 uv;
in float matIor;

uniform int   wire  = 0;
uniform float red   = 1.0;
uniform float green = 1.0;
uniform float blue  = 1.0;
uniform float alpha = 1.0;
uniform sampler2D sampler;
uniform vec2 widthBounds;
uniform vec2 lengthBounds;
//uniform float test = 0;

float rand(vec2 n) {
    return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}

float rand(float n) {
    return fract(sin(n) * 43758.5453123);
}

float rand(vec4 n) {
//    vec2 first2 = vec2(n[0], n[1]);
////    return 1.f;
    return rand(vec2(n[0] * rand(n[2]), n[1] * rand(n[3])));
}

vec2 uvFromWorldPoint(vec3 point) {
    float u = (point.x - widthBounds[0]) / (widthBounds[1] - widthBounds[0]);
    float v = (point.z - lengthBounds[0]) / (lengthBounds[1] - lengthBounds[0]);
    return vec2(u, v);
}

void main() {
    // Do lighting in camera space
    vec3 lightDir = -normalize(vec3(1, -1, 1));
//    lightDir = normalize(vec3(0.f, 3.f, 0.f) - pos);
//    float d = clamp(dot(normal_cameraSpace, lightDir), 0, 1);
    float d = clamp(dot(normal_worldSpace, lightDir), 0, 1);
    vec3 reflectedLight = lightDir - 2 * dot(lightDir, normal_worldSpace) * normal_worldSpace;
    vec3 posToCam = normalize(camera_worldSpace - pos);
    float spec = pow(clamp(dot(posToCam, reflectedLight), 0, 1), 2.f);

//    fragColor = texture(sampler, vec2(0.5f, 0.5f));
//    fragColor = vec4(abs(pos.x / 160.f), pos.y, 0.f, 1.f);
//    fragColor = vec4(uv.y, uv.y, 0.f, 1.f);
//    fragColor = vec4(camera_worldSpace.x - pos[0], camera_worldSpace.y - pos[1], pos[2], 1.f);
//    fragColor = vec4(- pos[0], 0.f, 0.f, 1.f);
//    fragColor = vec4((pos - vec3(widthBounds[0], 0, lengthBounds[0])) / 5.f, 1.f);
//    fragColor = vec4(fragColor.x, 0.f, fragColor.z, 1.f);
//    fragColor = vec4(test, test, test, 1.f);
    vec2 refrUV = uvFromWorldPoint(refrPos);
    float waterMurkiness = 0.002f; // TODO: Make uniform
    vec3 waterVolumeColor = vec3(red * 0.1f, green * 0.2f, blue * 0.2f);
    float murkDiffuse = 0.3f;
    float murkAmbient = 0.8f;
    float beerAtt = exp(-length((pos - refrPos)) * waterMurkiness);
    // EXPLANATION: WHEN THE WATER IS NOT PERFECTLY CLEAR, IT WILL HAVE STUFF SCATTERING LIGHT UNDERNEATH THE SURFACE
    // SOME OF IT WILL BE DIFFUSELY LIT, BUT THERE WILL BE SUBSURFACE SCATTERING, ESTIMATED BY THE AMBIENT TERM

    vec4 diffuse = vec4(red * d, green * d, blue * d, 1.0f);
    vec4 specular = vec4(1, 1, 1, 1) * pow(spec, 10.f);
//    vec4 transmissive = vec4(vec3(refrUV, 1.f - refrUV.y), 1.f);
    float waterBlurriness = 0.f;
    vec2 refrUVBlurry = (1 - beerAtt) * vec2(rand(refrUV), rand(vec4(pos, d))) * waterBlurriness + refrUV;
    vec4 transmissive = texture(sampler, vec2(refrUVBlurry));
    vec4 murk = (vec4(waterVolumeColor * d * murkDiffuse + waterVolumeColor * murkAmbient, 1.0f));

//    refrProb *= beerAtt;

    fragColor = 0.75f * diffuse; // Diffuse
    fragColor = vec4(0,0,0,0);
    fragColor += 1.f * specular; // Specular TODO: Pass multiplications as uniforms.
    fragColor = clamp(fragColor, 0.f, 1.f); // Clamp
    fragColor *=  ((1 - refrProb) / 1.f);

    vec4 volumetric = beerAtt * transmissive;
    volumetric += (1 - beerAtt) * murk;

    fragColor += refrProb * volumetric;
//    fragColor = transmissive * refrProb;
    fragColor = vec4(vec3(fragColor), 1.5f);
    // Dividing refrProb by 2 just for heuristic. Want more phong to show through.
//    fragColor = clamp(fragColor, 0.f, 1.f);
//    fragColor = vec4(refrProb, 0.f, 0.f, 1.f);

    // TODO: ACTUAL LIGHTING MODEL SHOULD BE SOMETHING LIKE
    // VELOCITY * DIFFUSE
    // (1 - refrProb) * SPECULAR
    // refrProb * (BEER * TRANSMISSIVE + (1 - beerAtt) * VOLUME (which is somewhat diffuse too?))
    // Transmissive shouldn't just get darker, but blurrier as beer attenuation lowers.
//    fragColor = texture(sampler, vec2(refrUV));
//    fragColor = vec4(normal_worldSpace[0], 0, normal_worldSpace[1], 1.f);
//    fragColor = diffuse;
}

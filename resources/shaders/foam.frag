#version 330 core

in vec2 constants;
in vec2 dir;
in vec2 tex;
in vec3 pos;
in vec3 norm;
in vec3 camera_worldSpace;


uniform float time;
uniform float phaseC; // phase constant

uniform sampler2D halftone_texture;
uniform sampler2D foam_texture;

uniform vec2 widthBounds;
uniform vec2 lengthBounds;

uniform vec4 sunColor = vec4(1.5f, .7, .39f, 1.f);
uniform vec3 lightDir = -normalize(vec3(1, 0, -1));

out vec4 fragColor;

float getSaturation(vec2 k, vec2 xzPos, float adjWaveLength, float phaseC){
    //k = normalize(k);
    float result = dot(k, xzPos) * 3.14f / adjWaveLength;
    result = result + phaseC*time*.5f;
    result = -tan(result + 1.57f);
    result = exp(result) / 20.f;

    return result;

}



void main() {
    float height = pos.y;
   float saturation = constants[0];//getSaturation(dir, vec2(pos.x, pos.z), 200.f, constants[0]);
   vec4 m_uv = texture(halftone_texture, tex*.6);
   float m_threshold = (m_uv.r + m_uv.g + m_uv.b) / 3;

   // final rgba color at x,z pos
   vec4 h = vec4(0,0,1,1);
   if (saturation > m_threshold) h = vec4(1,1,1, 1);

   // add fading effect to bubble popping
   vec4 g = clamp(saturation - m_threshold, 0, 1) * h;

   // apply foam texture
   vec4 foam = texture(foam_texture, tex + time*.0003);

    vec4 lightFromSun = vec4(1.f, 1.f, 1.f, 1.f);
    vec4 foamAmbient = vec4(.5f, .5f, .5f, 1.f);

    vec3 reflectedLight = lightDir - 2 * dot(lightDir, norm) * norm;
    vec3 posToCam = normalize(camera_worldSpace - pos);
    float spec = pow(clamp(dot(posToCam, reflectedLight), 0, 1), 2.f);

    float diffuseComp = dot(normalize(norm), lightDir);
    lightFromSun = vec4(vec3(diffuseComp * sunColor), 1.f);

   vec4 j = vec4(0,0,0,0);
   if (saturation > m_threshold) j = g * foam * (lightFromSun + foamAmbient + spec * .6f);


    fragColor = j;
    // fragColor = vec4(vec3(g), 1);
   //  fragColor = vec4(vec3(saturation), 1);
}

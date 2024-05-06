#version 330 core

in vec4 saturation_const;
in vec2 dir;
in vec2 tex;
in vec3 pos;


uniform float time;
uniform sampler2D halftone_texture;
uniform vec2 widthBounds;
uniform vec2 lengthBounds;

out vec4 fragColor;

float getSaturation(vec2 k, vec2 xzPos, float adjWaveLength, float phaseC){
    float result = dot(k, xzPos) * 3.14f / adjWaveLength;
    result = result + phaseC*time*.5f;
    result = -tan(result) + 1.57f;
    result = exp(result) / 4.f;

    return result;

}



void main() {
   //float saturation = getSaturation(saturation_const[0], saturation_const[1],saturation_const[2],saturation_const[3]);
   vec4 color = texture(halftone_texture, tex);

   fragColor = vec4(vec3(color), 1);
}

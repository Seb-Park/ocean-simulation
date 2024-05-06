#version 330 core

in vec2 constants;
in vec2 dir;
in vec2 tex;
in vec3 pos;


uniform float time;
uniform float phaseC; // phase constant

uniform sampler2D halftone_texture;
uniform sampler2D foam_texture;

uniform vec2 widthBounds;
uniform vec2 lengthBounds;

out vec4 fragColor;

float getSaturation(vec2 k, vec2 xzPos, float adjWaveLength, float phaseC){
    //k = normalize(k);
    float result = dot(k, xzPos) * 3.14f / adjWaveLength;
    result = result + phaseC*time*.5f;
    result = -tan(result + 1.57f);
    result = exp(result) / 4.f;

    return result;

}



void main() {
    float height = pos.y;
   float saturation = getSaturation(dir, vec2(pos.x,pos.z), 100.f, constants[0]);
   vec4 m_uv = texture(halftone_texture, tex);
   float m_threshold = m_uv.r * m_uv.g * m_uv.b;

   // final rgba color at x,z pos
   vec4 h = vec4(0,0,1,1);
   if (saturation > m_threshold) h = vec4(vec3(m_threshold), 1);

   // add fading effect to bubble popping
   vec4 g = clamp(.5*saturation - m_threshold, 0, 1) * h;

   // apply foam texture
   vec4 foam = texture(foam_texture, tex);
   vec4 j = vec4(0,0,0,1);
   if (saturation > m_threshold) j = vec4(vec3(g*foam), 1);


   fragColor = vec4(vec3(saturation), 1);
}

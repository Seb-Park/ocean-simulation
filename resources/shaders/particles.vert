#version 330 core
layout (location = 0) in vec4 pos; // <vec2 position, vec2 texCoords>

out vec2 TexCoords;
out vec4 ParticleColor;

uniform mat4 view, projection, model;
uniform vec3 offset;
uniform vec4 color;
uniform float alpha;


void main()
{
    float scale = 1.f;//2000.0f;
    TexCoords = vec2(pos.z, pos.w);
    //ParticleColor = color;
   // gl_Position = projection *view* vec4((pos * scale) + vec2(offset), 0.0, 1.0);
   gl_Position = (vec4(pos.x*alpha, pos.y*alpha, 0, 1) + projection*view*vec4(vec3(offset), 1))*scale;
   //gl_Position = vec4(pos.x, pos.y, 0, 1);


}

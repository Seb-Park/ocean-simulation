#version 330 core
in vec2 TexCoords;
in vec4 ParticleColor;

out vec4 fragColor;

uniform sampler2D sprite;
uniform float alpha;


void main()
{
   // color = (texture(sprite, TexCoords) * ParticleColor);
    fragColor = vec4(0,1,0,alpha);
}

#version 330 core
in vec2 TexCoords;
in vec4 ParticleColor;

out vec4 fragColor;

uniform sampler2D particle_texture;
uniform float alpha;


void main()
{
    vec4 c = texture(particle_texture, TexCoords);
    fragColor = c*vec4(1,1,1,alpha);
    //fragColor = vec4(1);
}

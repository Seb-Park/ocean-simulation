#version 330 core
out vec4 fragColor;
uniform sampler2D grass_texture;
uniform vec3 skyColor;


in GS_OUT {
    vec2 texCoord;
    float visibility;
} fs_in;

void main(void)
{
    vec4 color = texture(grass_texture, fs_in.texCoord);
    if (color.a < 0.4) discard;


    fragColor = color;
    fragColor = mix(vec4(skyColor, 1.f), fragColor, fs_in.visibility);
   // fragColor = vec4(0, 0, 1, 1);

}

#version 330 core
in vec2 texCoords;
out vec4 color;

uniform sampler2D guiTexture0;
uniform bool hovering;



void main()
{
    if (texCoords.x < 0.0 || texCoords.x > 1.0 ||
            texCoords.y < 0.0 || texCoords.y > 1.0) {
            discard;
    }

    vec4 texel0, texel1;
    //color = texture(guiTexture0, texCoords);

    texel0 = texture(guiTexture0, texCoords);
    //texel1 = texture(guiTexture1, texCoords);

    color = texel0;

    if (hovering){
        color = mix(color, vec4(0.f, 0.f, 0.f, 1.f), .5f);
    }

}

#version 330 core
in vec2 texCoords;
out vec4 color;

uniform sampler2D inventoryTexture0;

void main()
{
    if (texCoords.x < 0.0 || texCoords.x > 1.0 ||
            texCoords.y < 0.0 || texCoords.y > 1.0) {
            discard;
    }

    color = texture(inventoryTexture0, texCoords);
}

#version 330 core
layout (location = 0) in vec2 pos; // <vec2 pos, vec2 tex>
out vec2 texCoords;
uniform mat4 transform;
uniform mat4 projection;
uniform vec2 textureScale;
uniform bool hovering;




void main()
{

    gl_Position = transform*vec4(pos.x, pos.y, 0.f, 1.f);
    texCoords = textureScale*vec2((pos.x+1.f)/2.f, (pos.y+1.f)/2.f);
    //texCoords =  textureScale * (pos.x - 0.5f) + 0.5f;
}

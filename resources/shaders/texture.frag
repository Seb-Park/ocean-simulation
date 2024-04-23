#version 330 core

// TASK 16: Create a UV coordinate in variable
in vec2 uv;

// TASK 8: Add a sampler2D uniform
uniform sampler2D sampler;

// TASK 29: Add a bool on whether or not to filter the texture
uniform bool filtered;

out vec4 fragColor;

void main()
{
    // TASK 17: Set fragColor using the sampler2D at the UV coordinate
    fragColor = texture(sampler, uv);

    // TASK 33: Invert fragColor's r, g, and b color channels if your bool is true
    if(filtered){
        fragColor = vec4(1) - fragColor;
        fragColor.w = 1;
    }
}

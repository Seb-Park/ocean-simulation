#version 330 core
out vec4 fragColor;

// Additional information for lighting
in vec4 normal_worldSpace;
in vec4 position_worldSpace;
in vec4 force_worldSpace;

uniform int wire = 0;
uniform float red = 1.0;
uniform float green = 1.0;
uniform float blue = 1.0;
uniform float alpha = 1.0;
uniform int displayForce = 0;

void main() {
    if (wire == 1) {
        fragColor = vec4(0.0, 0.0, 0.0, 1);
        return;
    }
    vec4 lightPos   = vec4(-2.0, 2.0, -3.0 , 1.0);
    vec3 lightColor = vec3(1.0f, alpha, 0.0f);
    vec4 lightDir   = normalize(-lightPos + position_worldSpace);
    float c = clamp(dot(-normal_worldSpace, lightDir), 0, 1);

    if(displayForce == 1) {
        fragColor = force_worldSpace;
        return;
    }

    float r = red;
    float g = green;
    float b = blue;

//    if (displayForce == 1 && wire == 0) {
//        r = force_worldSpace[0];
//        g = force_worldSpace[1];
//        b = force_worldSpace[2];
//    }

    fragColor = vec4(r * c * lightColor[0], g * c * lightColor[0], b * c * lightColor[0], 1);

}

#version 330 core
// Uniforms for shape information
in vec3 worldSpace_pos;
in vec3 worldSpace_norm;
in vec2 tex_coord;

in float visibility;

uniform vec3 skyColor;

// Object Material Data
uniform int colorSource; // 0 = solid color (objColor), 1 = texture color (objTexture), 2 = per-vertex color (vertColor)
uniform vec3 objColor;
uniform sampler2D objTexture;
in vec3 vertColor;
uniform float shininess;

// Camera uniform
uniform vec3 worldSpace_camPos;

// Global Data
uniform vec3 coeffs; // vec3(ka, kd, ks)

// Light Data
uniform int lightType[16]; // 0 = point light, 1 = directional light
uniform vec3 lightColor[16];
uniform vec3 lightFunction[16]; // Attenuation coefficients
uniform vec3 worldSpace_lightPos[16]; //Light Positions
uniform vec3 worldSpace_lightDir[16]; //Light Directions
uniform int numLights; // Max number of lights = 8

out vec4 fragColor;

vec3 getToLight(int lightIndex) {
    int LIGHT_POINT = 0;
    int LIGHT_DIRECTIONAL = 1;

    if (lightType[lightIndex] == LIGHT_POINT) {
        return normalize(worldSpace_lightPos[lightIndex] - worldSpace_pos);
    }
    else if (lightType[lightIndex] == LIGHT_DIRECTIONAL) {
        return normalize(-worldSpace_lightDir[lightIndex]);
    }

    return vec3(0);
}

float attenuationFactor(int lightIndex) {
    int LIGHT_POINT = 0;

    if (lightType[lightIndex] == LIGHT_POINT) {
        vec3 coeffs = lightFunction[lightIndex];
        float d = length(worldSpace_lightPos[lightIndex] - worldSpace_pos);
        return 1.0 / (coeffs.x + coeffs.y * d + coeffs.z * d * d);
    }

    return 1;
}

float computeDiffuseIntensity(vec3 worldSpace_toLight) {
    // Dot product to get diffuse intensity
    return max(dot(worldSpace_toLight, normalize(worldSpace_norm)), 0);
}

float computeSpecularIntensity(vec3 worldSpace_toLight, vec3 worldSpace_toEye) {
    // Guard against pow weirdness when exponent is 0
    if (shininess == 0) {
        return 0;
    }

    //reflect toLight
    vec3 worldSpace_toLightReflected = reflect(-worldSpace_toLight, normalize(worldSpace_norm));

    //Compute specular intensity using toEye, reflected light, and shininess
    return pow(max(dot(worldSpace_toLightReflected, worldSpace_toEye), 0), shininess);
}

void main() {
    // Declare ambient, diffuse, and specular terms
    vec3 ambi = vec3(coeffs.x);
    vec3 diff = vec3(0.0);
    vec3 spec = vec3(0.0);


    // Compute worldSpace_toEye Vector for specular intensity computation;
    vec3 worldSpace_toEye = normalize(worldSpace_camPos - worldSpace_pos);


    // Compute per-light diffuse and specular contribution
    for(int i = 0; i<numLights; i+= 1){

        // get direction vector to light based on light type
        vec3 worldSpace_toLight = getToLight(i);

        float diffuse_intensity = computeDiffuseIntensity(worldSpace_toLight);
        float specular_intensity = computeSpecularIntensity(worldSpace_toLight, worldSpace_toEye);

        float att = attenuationFactor(i);


        diff = diff + diffuse_intensity * lightColor[i] * att;
        spec = spec + specular_intensity * lightColor[i] * att;
    }

    // Apply global coefficients and object color to the diffuse and specular components
    diff = diff * vec3(coeffs.y);
    spec = spec * vec3(coeffs.z);

    // Color generated only from light intensities and colors
    vec3 tempColor = clamp(ambi + diff + spec, 0, 1);

    // Apply correct object color
    if (colorSource == 0 ) {
        fragColor = vec4(tempColor * objColor, 1.0);
    } 
    else if (colorSource == 1){
        fragColor = vec4(tempColor * vec3(texture(objTexture, tex_coord)), 1.0);
    }
    else if (colorSource == 2) {
        fragColor = vec4(tempColor * vertColor, 1.0);
    }
    else{
        fragColor = vec4(tempColor, 1.0);
    }

    // mix actual color with fog color
   // fragColor = mix(vec4(.77f, .85f, .99f, 1.f), fragColor, visibility);
}

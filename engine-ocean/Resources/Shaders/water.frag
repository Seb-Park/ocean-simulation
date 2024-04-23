#version 330 core
// Uniforms for shape information
//in vec3 worldSpace_pos;
//in vec3 worldSpace_norm;

in float visibility;
in vec4 clipSpace;
in vec2 tex_coord;
in vec3 toCameraVector;



// Camera uniform
uniform vec3 worldSpace_camPos;
uniform sampler2D reflect_texture;
uniform sampler2D refract_texture;
uniform sampler2D du_dv_map;
uniform float moveFactor; // offset that changes over time

out vec4 fragColor;

const float distortionStrength = .02f;
const vec4 waterColor = vec4(.06f, .49f, .59f, 1.f);
uniform vec3 skyColor;


void main() {

    // screen space coords
    vec2 ndc = (clipSpace.xy/clipSpace.w)/2.f + .5f;

    vec2 refractTexCoords = ndc;
    vec2 reflectTexCoords = vec2(ndc.x, ndc.y);

    // get associated rg color in the range -1,1 from du_dv map
    vec2 distortion1 = (texture(du_dv_map, vec2(tex_coord.x + moveFactor, tex_coord.y)).rg*2.f - 1.f) * distortionStrength;
    vec2 distortion2 = (texture(du_dv_map, vec2(-tex_coord.x + moveFactor, tex_coord.y + moveFactor)).rg*2.f - 1.f) * distortionStrength;
    vec2 totalDistortion = distortion1 + distortion2;

    refractTexCoords += totalDistortion;
    reflectTexCoords += totalDistortion;

    refractTexCoords = clamp(refractTexCoords, 0.001f, .999f);
    reflectTexCoords = clamp(reflectTexCoords, 0.001f, .999f);

    // fresnel effect
    vec3 viewVector = normalize(toCameraVector);
    float refractiveFactor = dot(viewVector, vec3(0.f, 1.f, 0.f));
    refractiveFactor = pow(refractiveFactor, 5.f); // changes how reflective water is





    vec4 reflect = texture(reflect_texture, reflectTexCoords);
    vec4 refract = texture(refract_texture, refractTexCoords);


    // mix actual color with fog color
    fragColor = mix(reflect, refract, refractiveFactor);
    fragColor = mix(fragColor, waterColor, .2f);
    fragColor = mix(vec4(skyColor, 1.f), fragColor, visibility);
}

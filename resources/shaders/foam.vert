#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in float wavelength; // wavelenth adjusted for ocean depth
layout(location = 2) in vec2 direction; // wave slope
layout(location = 3) in vec2 texCoords; // texture coords
layout(location = 3) in vec3 norm; // texture coords


out vec4 saturation_const;
out vec2 dir;
out vec2 tex;
out vec3 pos;



uniform float time;
uniform float phaseC; // phase constant
uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;
uniform mat4 inverseView;
uniform vec2 widthBounds;
uniform vec2 lengthBounds;

vec2 calculateTexCoord(vec3 pos){
//    vec3 v = vec3(0,0,1);
//    if (abs(norm.y) < 1.f){
//        v = normalize(vec3(0,1,0) - norm.y*norm);
//    }
//    vec3 u = normalize(cross(norm, v));

//    float u_coord = dot(u, vec3(pos.x, 0, pos.z)) - widthBounds[0]/ (widthBounds[1] - widthBounds[0]);
//    float v_coord = dot(v, vec3(pos.x, 0, pos.z)) - lengthBounds[0]/ (lengthBounds[1] - lengthBounds[0]);

     float u_coord = position.x / (widthBounds[1] - widthBounds[0]);
     float v_coord = position.z / (lengthBounds[1] - lengthBounds[0]);


     float offset = .5f;
    return vec2(u_coord + offset, v_coord + offset);

}

void main() {
    dir = direction;
    saturation_const = vec4(wavelength, phaseC, position.x, position.z);

    gl_Position = proj * view * model * vec4(position, 1);
    pos = vec3(gl_Position);

    tex = calculateTexCoord(position);

}

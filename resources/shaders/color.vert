#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex
layout(location = 2) in vec3 texCoords;   // Normal of the vertex

uniform float depth = -1000.f;
uniform float skyHeight = 500.f;
uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;
uniform mat4 inverseView;
//uniform float width;
//uniform float height; // TODO: Pass in width and height as uniform

uniform mat3 inverseTransposeModel;

out vec3 normal_cameraSpace;
out vec3 normal_worldSpace;
out vec3 camera_worldSpace;
out vec3 pos;
out float intensity;
out vec3 oldPosFlat;

uniform vec2 widthBounds;
uniform vec2 lengthBounds;

vec3 moveToTopDown(vec3 p) {
    return vec3((p[0] - ((widthBounds[0] + widthBounds[1]) / 2)) / (widthBounds[1] - widthBounds[0]) * 2,
            -(p[2] - ((lengthBounds[0] + lengthBounds[1]) / 2)) * 2 / (lengthBounds[1] - lengthBounds[0]), 0.f);
}

vec4 refractToFloor(vec3 l, vec3 p, vec3 n, float d) {
    // Refracts incoming light direction l through normal n at point p until hits floor at depth d
    vec3 w_o = normalize(l);
    float cos_theta_i = dot(-w_o, n);
    float n_i = 1;
    float n_t = 1.33f;
    float determinant = 1.f - (pow((n_i / n_t), 2.f) * (1.f - pow(cos_theta_i, 2.f)));

    float r0 = pow((n_i - n_t) / (n_i + n_t), 2.f); // variable required to calculate probability of reflection
    float prob_to_refl = r0 + ((1 - r0) * pow((1 - cos_theta_i), 5.f));

    if (determinant >= 0) {
        float cos_theta_t = sqrt(determinant);
        vec3 w_t = (n_i / n_t) * w_o + ((n_i / n_t) * cos_theta_i - cos_theta_t) * n;
        float dist = p.y - d;
        float depthScale = dist / w_t.y;
        vec3 groundContactPoint = -(w_t * depthScale) + p;
        return vec4(groundContactPoint, 1.f - prob_to_refl);
    } else {
        return vec4(0, 0, 0, 0);
    }
}

vec3 wrapGL(vec3 p) {
    float newX = p[0];
    float newY = p[1];

    newX = mod(newX + 1.f, 2.f) - 1.f;

    newY = mod(newY + 1.f, 2.f) - 1.f;

    return vec3(newX, newY, 0.f);
}

void main() {
    normal_cameraSpace = normalize(inverse(transpose(mat3(view))) * inverseTransposeModel * normal);
    camera_worldSpace = vec3(inverseView * vec4(0.f, 0.f, 0.f, 1.f));
    normal_worldSpace = normal;
    pos = vec3(model * vec4(position, 1.f)); //vec3(model * vec4(objSpacePos, 1.f));
    vec3 lightDir = -normalize(vec3(1, -1, 1));

    gl_Position = proj * view * model * vec4(position, 1);

//    gl_Position = vec4((position[0] - ((widthBounds[0] + widthBounds[1]) / 2)) / (widthBounds[1] - widthBounds[0]) * 2,
//            (position[2] - ((lengthBounds[0] + lengthBounds[1]) / 2)) / (lengthBounds[1] - lengthBounds[0]) * 2, 0.f, 1.f);
    float newX = mod(int(position[0]), 3) - 1;
    float newY = mod(int(position[2]), 3) - 1;

    vec4 refractedPositionAndProb = refractToFloor(-lightDir, position, normal_worldSpace, depth);

    float waterMurkiness = .002f;

    gl_Position = vec4(newX, newY, 0.f, 1.f);

    intensity = refractedPositionAndProb[3] * clamp(dot(normal_worldSpace, lightDir), 0.f, 1.f) * exp(-length((position - vec3(refractedPositionAndProb))) * waterMurkiness);

    oldPosFlat = moveToTopDown(position);
    pos = moveToTopDown(vec3(refractedPositionAndProb));
    gl_Position = vec4(pos, 1.f);
    gl_Position = vec4(oldPosFlat, 1.f);
}

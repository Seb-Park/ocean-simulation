#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex
layout(location = 2) in vec3 texCoords;   // Normal of the vertex

uniform float depth = -1000.f;
uniform float skyHeight = 500.f;

uniform sampler2D normSamp;

out vec3 normal_worldSpace;
out vec3 pos;
out vec3 newPos;
out vec4 col;
out float refractProb;

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
        float dist = p.z - d;
        float depthScale = dist / w_t.z;
//        vec3 groundContactPoint = -(w_t * depthScale) + p;
        vec3 groundContactPoint = (w_t * depthScale) + p;
        return vec4(groundContactPoint, 1.f - prob_to_refl);
    } else {
        return vec4(0, 0, 0, 0);
    }
}

void main() {
    normal_worldSpace = normal;
    pos = position;
    vec4 sampledNormal = texture(normSamp, vec2((pos + 1) / 2));
    sampledNormal = (sampledNormal * 2.f) - 1.f;
    col = sampledNormal;
    vec4 newPosAndProb = refractToFloor(vec3(0, 0, 1), pos, normalize(vec3(sampledNormal)), 0.01f);
    newPos = vec3(newPosAndProb[0], newPosAndProb[1], 0.f);
    refractProb = newPosAndProb[3];
//    newPos = pos;
    gl_Position = vec4(newPos, 1.f);
}

#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex
layout(location = 3) in vec3 texCoords;   // Normal of the vertex


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
out vec3 reflPos;
out vec3 refrPos;
out float refrProb;
out vec2 uv;
out float matIor;

vec4 getRefrPos() {
    float depth = -1000.f; // TODO: Pass as uniform
    vec3 w_o = normalize(pos - camera_worldSpace);
    float cos_theta_i = dot(-w_o, normal_worldSpace);
    float n_i = 1;
    float n_t = 1.33f;
//    matIor = n_t;
    float determinant = 1.f - (pow((n_i / n_t), 2.f) * (1.f - pow(cos_theta_i, 2.f)));

    float r0 = pow((n_i - n_t) / (n_i + n_t), 2.f); // variable required to calculate probability of reflection
    float prob_to_refl = r0 + ((1 - r0) * pow((1 - cos_theta_i), 5.f));

    if (determinant >= 0) {
        float cos_theta_t = sqrt(determinant);
        vec3 w_t = (n_i / n_t) * w_o + ((n_i / n_t) * cos_theta_i - cos_theta_t) * normal_worldSpace;
//        Ray reflectedRay(i.hit, w_t);
//        float attenuation = (!entering && attenuateRefract) ? std::pow(M_E, (-(ray.o - i.hit).norm()) * mat.ior) : 1;
//        L += traceRay(reflectedRay, scene, true, n_t, !is_in_refractor) * attenuation / threshold;
        float dist = position.y - depth;
        float depthScale = dist / w_t.y;
        vec3 groundContactPoint = -(w_t * depthScale) + position;
        return vec4(groundContactPoint, 1.f - prob_to_refl);
    } else {
//        Eigen::Vector3f w_i = w_o - 2 * w_o.dot(intersectNormal) * incidenceNormal;
//        Ray reflectedRay(i.hit, w_i);
//        L += traceRay(reflectedRay, scene, true, current_ior, false) / threshold;
        return vec4(0, 0, 0, 0);
    }
}

vec3 getReflPos() {
    float depth = 5000.f; // TODO: Pass as uniform
    vec3 w_o = normalize(pos - camera_worldSpace);
    vec3 reflectedRay = 2 * dot(-w_o, normal_worldSpace) * normal_worldSpace + w_o;

    float dist = depth - position.y;
    float depthScale = dist / reflectedRay.y;
    vec3 skyContactPoint = (reflectedRay * depthScale) + position;
    return skyContactPoint;
}

void main() {
//    float depth = -4.f;
//    float dist = position.y - depth;
    float width = 81.f * 2.f;
    float length = 81.f * 2.f;
    matIor = 1.33f;

    normal_cameraSpace = normalize(inverse(transpose(mat3(view))) * inverseTransposeModel * normal);
    camera_worldSpace = vec3(inverseView * vec4(0.f, 0.f, 0.f, 1.f));
    normal_worldSpace = normal;
    pos = vec3(model * vec4(position, 1.f)); //vec3(model * vec4(objSpacePos, 1.f));
//    pos = position;

//    float depthScale = dist / normal.y;
//    vec3 groundContactPoint = -(normal * depthScale) + position; // carries down to ground
//    groundContactPoint = vec3(model * vec4(position, 1));
//    uv = vec2((position.x + 81.f) / (162.f), groundContactPoint.z);
//    uv = vec2(normal);
    vec4 refrPos_and_prob = getRefrPos();
    refrPos = vec3(refrPos_and_prob);
    refrProb = clamp(refrPos_and_prob.w, 0.f, 1.f);

    reflPos = getReflPos();

    gl_Position = proj * view * model * vec4(position, 1);
}

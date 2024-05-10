#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include <iostream>
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE
#include "graphics/shader.h"
#include "graphics/camera.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

struct Particle{
    Eigen::Vector3f pos = Eigen::Vector3f(0,0,0);
    Eigen::Vector3f vel = Eigen::Vector3f(0,2000,0);

    Eigen::Vector4f color = Eigen::Vector4f(1,1,1,1);
    float life = 1.f;

};

class particlesystem
{
public:
    particlesystem();

    void update(double deltaTime);
    void draw(Shader *skybox_shader, Camera  m_camera);
    void draw(Shader *shader, Camera  m_camera, std::vector<Eigen::Vector3f> verts, Eigen::Matrix4f model);

    void init();

    void setVerts(std::vector<Eigen::Vector3f> verts){
        std::cout << verts.size() << std::endl;
        m_verts = verts;

    }



private:
    int m_amount = 500;
    int m_new_amount = 100; // new particles per frame
    int lastUsedIndex = 0;
    Eigen::Vector3f gravity = Eigen::Vector3f(0,-9.81f,0);
    float factor = 100.f; // velocity multiplier



    void respawn_particle(Particle &p, Eigen::Vector3f new_pos);
    int getUnusedParticleIndex();
    Eigen::Vector3f getRandomInitialVel();




    float d = 100.f;
    std::vector<float> m_vertices = {
        -d, d,
        d,d,
        -d, -d,
        d, d,
        d, -d,
        -d, -d
    };

    GLuint VAO, VBO; // holds quad shape






    std::vector<Particle> m_particles;

    std::vector<Eigen::Vector3f> m_verts;

};

#endif // PARTICLESYSTEM_H

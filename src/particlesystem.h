
#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include "ocean/ocean_alt.h"
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
    void draw(Shader *shader, Camera  m_camera, GLuint texture);

    void init(std::vector<OceanSpray> verts);

    void setVerts(std::vector<OceanSpray> verts){
        //std::cout << "VERTS SIZE:" << verts.size() << std::endl;
        m_verts = verts;

    }



private:
    int m_amount = 500;
    int m_new_amount = 100; // new particles per frame
    int lastUsedIndex = 0;
    Eigen::Vector3f gravity = Eigen::Vector3f(0,-9.81f,0);
    float factor = 40.f; // velocity multiplier



    void respawn_particle(Particle &p, OceanSpray new_pos);
    int getUnusedParticleIndex();
    Eigen::Vector3f getRandomInitialVel();




    float d = 50.f;
    std::vector<float> m_vertices = {
        -d, d, 0,1,
        -d, -d, 0, 0,
        d,d, 1,1,
        //
        d,d, 1,1,
        -d, -d, 0, 0,

        d, -d, 1, 0
    };

    GLuint VAO, VBO; // holds quad shape






    std::vector<Particle> m_particles;

    std::vector<OceanSpray> m_verts;

};

#endif // PARTICLESYSTEM_H

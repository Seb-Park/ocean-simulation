#include "particlesystem.h"
#include <iostream>

particlesystem::particlesystem()
{

}

std::pair<double,double> sample_complex_gaussian(){
    double uniform_1 = (double)rand() / (RAND_MAX);
    double uniform_2 = (double)rand() / (RAND_MAX);

    // set a lower bound on zero to avoid undefined log(0)
    if (uniform_1 == 0)
    {
        uniform_1 = 1e-10;
    }
    if (uniform_2 == 0)
    {
        uniform_2 = 1e-10;
    }

    // real and imaginary parts of the complex number
    double real = sqrt(-2*log(uniform_1)) * cos(2*M_PI*uniform_2);
    double imag = sqrt(-2*log(uniform_1)) * sin(2*M_PI*uniform_2);

    return std::make_pair(real, imag);
}

float getRandomInRange(float max, float min){
    return min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));

}

float getRandomY(){
    return getRandomInRange(.001,.2);
}
Eigen::Vector3f particlesystem::getRandomInitialVel(){


    std::pair<float, float> r = sample_complex_gaussian();
    return factor*Eigen::Vector3f(r.first, getRandomY(), r.second);
}

void particlesystem::init(std::vector<OceanSpray> verts){
    // make sure to set m_verts
    setVerts(verts);

    m_particles.reserve(4000);

    for (auto v : m_verts){
        Particle p;
        p.pos = Eigen::Vector3f(v.height);
        p.vel = v.slope*factor;
        //p.vel = getRandomInitialVel();
        p.vel = getRandomInitialVel().asDiagonal() * v.slope;
        p.vel[1] = 0.1;


        m_particles.push_back(p);
    }


    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, m_vertices.size()*sizeof(float), m_vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0);

    // unbind
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void particlesystem::update(double deltaTime){

    // add new particles
    for (int i=0; i<m_particles.size(); i++){
        int particle_index = getUnusedParticleIndex();
        if (particle_index == -1) continue;
        if (particle_index >= m_particles.size() - 1 || particle_index >= m_verts.size() - 1) continue;
        respawn_particle(m_particles[particle_index], m_verts[particle_index]);

    }


    float dt = deltaTime;
    // update all particles values
    for (Particle &p : m_particles){
       p.life -= deltaTime * getRandomInRange(.1f, .5f);

        // if particle is still alive, update pos
        if (p.life >= 0.f){
            float r = getRandomInRange(.5, 1.f);
            p.vel += gravity*dt*r;
            p.pos += p.vel * dt * 10.f;
//            p.vel[1] -= (p.vel.dot(p.vel) / 2 + gravity[0]) * dt;
//            p.pos += p.vel * dt * 10.f;
        }

        if (p.pos[1] < 0.f) p.life = 0.f;

    }
}



// finds the first instance of a particle that is dead, and returns its index so it can be replaced
int particlesystem::getUnusedParticleIndex(){

    // search from last used index
    for (int i=lastUsedIndex; i<m_verts.size(); ++i){
        if (m_particles[i].life <= 0.f){
            lastUsedIndex = i;
            return i;
        }
    }

    // otherwise do linear search up to last used
    for (int i=0; i<lastUsedIndex; ++i){
        if (m_particles[i].life <= 0.f){
            lastUsedIndex = i;
            return i;
        }
    }

    // if no dead ones found, don't respawn
    lastUsedIndex = 0;
    //std::cout << "first particle override!" << std::endl;
    return -1;
}

void particlesystem::respawn_particle(Particle &p, OceanSpray new_pos){

    p.pos = Eigen::Vector3f(new_pos.height);
    p.vel = new_pos.slope*factor;
    //p.vel = getRandomInitialVel();
    p.vel = getRandomInitialVel().asDiagonal() * new_pos.slope;
   p.vel[1] = .01;


    // reset life
    p.life = 1.f;

}

void particlesystem::draw(Shader *shader, Camera  m_camera){

    shader->bind();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE);


    // activate texture
//    glActiveTexture(GL_TEXTURE10);
//    glBindTexture(GL_TEXTURE_CUBE_MAP, skybox_tex);
//    glUniform1i(glGetUniformLocation(shader->id(), "cubeMap"), 9); // bind texture at slot 9


    // manually set view and projection, for non-translating view
    Eigen::Matrix4f projection = m_camera.getProjection();
    Eigen::Matrix4f view = m_camera.getView();
   // view.col(3) = Eigen::Vector4f(0, 0, 0, 0);

    glUniformMatrix4fv(glGetUniformLocation(shader->id(), "view"), 1, GL_FALSE, view.data());
    glUniformMatrix4fv(glGetUniformLocation(shader->id(), "projection"), 1, GL_FALSE, projection.data());
   // shader->setUniform("model", model);

    int i = 0;
    for (Particle p : m_particles){
        if (p.life >= 0.f){
            shader->setUniform("offset", p.pos);
            shader->setUniform("alpha", p.life);

            glBindVertexArray(VAO);
            glDrawArrays(GL_TRIANGLES, 0, m_vertices.size());
            glBindVertexArray(0);
            i++;
        }
    }
    shader->unbind();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


}

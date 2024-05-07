#include "skybox.h"

#include "stb/stb_image.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>





skybox::skybox()
{
    //initializeVAO();
}

void skybox::initializeVAO(){
    sky_shape.init(m_vertices_eigen, m_faces);
    std::cout << "hehee" << std::endl;
    skybox_tex = loadCubeMap(m_skyboxTextureFiles);

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, m_vertices.size()*sizeof(float), m_vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);

    // unbind
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void skybox::draw(Shader *skybox_shader, Camera  m_camera){




    glDepthFunc(GL_LEQUAL);
    glDisable(GL_CULL_FACE);

    skybox_shader->bind();
    Eigen::Vector3f sc = Eigen::Vector3f(.77f, .85f, .99f); // skycolor for fade effect
    glUniform3f(glGetUniformLocation(skybox_shader->id(), "skyColor"), sc[0], sc[1], sc[2]);


    // activate texture
    glActiveTexture(GL_TEXTURE9);
    glBindTexture(GL_TEXTURE_CUBE_MAP, skybox_tex);
    glUniform1i(glGetUniformLocation(skybox_shader->id(), "cubeMap"), 9); // bind texture at slot 9


    // manually set view and projection, for non-translating view
    Eigen::Matrix4f projection = m_camera.getProjection();
    Eigen::Matrix4f view = m_camera.getView();
    view.col(3) = Eigen::Vector4f(0, 0, 0, 0);

    glUniformMatrix4fv(glGetUniformLocation(skybox_shader->id(), "view"), 1, GL_FALSE, view.data());
    glUniformMatrix4fv(glGetUniformLocation(skybox_shader->id(), "projection"), 1, GL_FALSE, projection.data());

//    skybox_shader->setUniform("projection", projection);
//    skybox_shader->setUniform("view", view);

    // apply rotation matrix
    skybox_shader->setUniform("rotation", m_rotation_mat);

    //glUniformMatrix4fv(glGetUniformLocation(skybox_shader., "rotation"), 1, GL_FALSE, glm::value_ptr(m_rotation_mat->getRotation()[0]));

    //Global::graphics.setGlobalData(glm::vec3(.5f));
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, m_vertices.size());


//    glDrawElements(GL_TRIANGLES, m_vertices.size(), GL_UNSIGNED_INT, reinterpret_cast<GLvoid *>(0));
//    glBindVertexArray(0);
////        glBindTexture(GL_TEXTURE_2D, 0);
///
///
///
///

    std::cout << m_camera.getPosition() << std::endl;

   // sky_shape.draw(skybox_shader, GL_TRIANGLES);

   glDepthFunc(GL_LESS);


          skybox_shader->unbind();

}

void skybox::update(double deltaTime){
    m_rotation += ROTATE_SPEED * deltaTime;

//    Eigen::Rotation2Df f = Eigen::Rotation2Df(m_rotation);


//    m_rotation_mat = f.toRotationMatrix();


        auto sinA = std::sin(m_rotation / 2);
        auto cosA = std::cos(m_rotation/ 2);

        Eigen::Quaternionf q;
        q.x() = 0 * sinA;
        q.y() = 1 * sinA;
        q.z() = 0 * sinA;
        q.w() = cosA;

        Eigen::Matrix3f mat3 = q.toRotationMatrix();
        Eigen::Matrix4f mat4 = Eigen::Matrix4f::Identity();
        mat4.block(0,0,3,3) = mat3;


        m_rotation_mat = mat4;






}

GLuint skybox::loadCubeMap(std::vector<const char*> textureFiles){
    // create empty texture
    GLuint textureID;
    glGenTextures(1, &textureID);

    std::cout << "hello fssd" << std::endl;


    //glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    std::cout << "hello fssd" << std::endl;

    GLuint target = GL_TEXTURE_CUBE_MAP_POSITIVE_X;
    for (int i=0; i<6; i++){
        std::string filename = std::string(textureFiles[i]);//directory + '/' + filename;
        int width, height, nrChannels;
        unsigned char* data = stbi_load(filename.c_str(), &width, &height, &nrChannels, 0);

        if (data){
            stbi_set_flip_vertically_on_load(false);
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);

        }    else {
            std::cout << "Texture failed to load at path: " << textureFiles[i] << std::endl;
            stbi_image_free(data);
        }
    }

    return textureID;
}

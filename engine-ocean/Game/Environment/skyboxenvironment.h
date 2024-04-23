#ifndef SKYBOXENVIRONMENT_H
#define SKYBOXENVIRONMENT_H
#include "Graphics/global.h"
#include "Graphics/camera.h"
#include <GLFW/glfw3.h>
#include "Environment.h"


class SkyboxEnvironment : public Environment
{
public:
    SkyboxEnvironment(std::shared_ptr<Camera> camera);
    void draw() override;
    void update(double deltaTime) override;
private:
    void initializeVAO();

    GLuint VAO;
    GLuint VBO;
    GLuint EBO;
    GLuint skybox_tex;
    std::shared_ptr<Camera> m_camera;

    float SIZE = 500.f;

    std::vector<int> m_vertex_indices = {
        // Right
            1, 2, 6,
            6, 5, 1,
            // Left
            0, 4, 7,
            7, 3, 0,
            // Top
            4, 5, 6,
            6, 7, 4,
            // Bottom
            0, 3, 2,
            2, 1, 0,
            // Back
            0, 1, 5,
            5, 4, 0,
            // Front
            3, 7, 6,
            6, 2, 3
    };
    std::vector<GLfloat> m_vertices = {
        -SIZE,  SIZE, -SIZE,
        -SIZE, -SIZE, -SIZE,
         SIZE, -SIZE, -SIZE,
         SIZE, -SIZE, -SIZE,
         SIZE,  SIZE, -SIZE,
        -SIZE,  SIZE, -SIZE,

        -SIZE, -SIZE,  SIZE,
        -SIZE, -SIZE, -SIZE,
        -SIZE,  SIZE, -SIZE,
        -SIZE,  SIZE, -SIZE,
        -SIZE,  SIZE,  SIZE,
        -SIZE, -SIZE,  SIZE,

         SIZE, -SIZE, -SIZE,
         SIZE, -SIZE,  SIZE,
         SIZE,  SIZE,  SIZE,
         SIZE,  SIZE,  SIZE,
         SIZE,  SIZE, -SIZE,
         SIZE, -SIZE, -SIZE,

        -SIZE, -SIZE,  SIZE,
        -SIZE,  SIZE,  SIZE,
         SIZE,  SIZE,  SIZE,
         SIZE,  SIZE,  SIZE,
         SIZE, -SIZE,  SIZE,
        -SIZE, -SIZE,  SIZE,

        -SIZE,  SIZE, -SIZE,
         SIZE,  SIZE, -SIZE,
         SIZE,  SIZE,  SIZE,
         SIZE,  SIZE,  SIZE,
        -SIZE,  SIZE,  SIZE,
        -SIZE,  SIZE, -SIZE,

        -SIZE, -SIZE, -SIZE,
        -SIZE, -SIZE,  SIZE,
         SIZE, -SIZE, -SIZE,
         SIZE, -SIZE, -SIZE,
        -SIZE, -SIZE,  SIZE,
         SIZE, -SIZE,  SIZE
    };

    std::vector<const char*> m_skyboxTextureFiles = {"/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/skybox/hills2/right.png",
                                                  "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/skybox/hills2/left.png",
                                                  "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/skybox/hills2/up.png",
                                                  "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/skybox/hills2/bottom.png",
                                                  "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/skybox/hills2/back.png",
                                                  "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/skybox/hills2/front.png"};

    float ROTATE_SPEED = .1f; // 1 degree per sec
    float m_rotation = 0.f;
    std::shared_ptr<ModelTransform> m_rotation_mat;

};

#endif // SKYBOXENVIRONMENT_H

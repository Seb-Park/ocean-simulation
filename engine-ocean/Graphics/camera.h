#pragma once

#include "glm/glm.hpp"
#include "glm/ext.hpp"

class Camera
{
public:
    Camera(int width = 640, int height = 480,
    glm::vec3 pos = glm::vec3(0, 0, 0), glm::vec3 look = glm::vec3(0, 0, 1), 
    glm::vec3 up = glm::vec3(0, 1, 0), float fov = 1.f, 
    float nearPlane = 0.1f, float farPlane = 100.f);
    ~Camera();

    // Functions to get camera data for drawing
    glm::mat4 getProjection();
    glm::mat4 getView();

    // Functions to edit camera
    void resize(int width, int height);
    void translate(glm::vec3 move);
    void setPos(glm::vec3 newPos);
    glm::vec3 getPos();
    void rotate(float angle, glm::vec3 axis);
    void setLook(glm::vec3 newLook);
    glm::vec3 getLook();
    void setUp(glm::vec3 newUp);
    glm::vec3 getUp();

    glm::vec3 getLookAt();
    glm::vec3 getViewDirection();
    glm::vec3 getRight();


    int getHeight();
    int getWidth();



private:
    // Internal functions to update projection and view matrices
    void calculateProjection();
    void calculateView();

    int m_width;
    int m_height;
    glm::vec3 m_pos;
    glm::vec3 m_look;
    glm::vec3 m_up;
    float m_fov;
    float m_aspect;
    float m_near;
    float m_far;

    glm::mat4 m_proj = glm::mat4(1);
    glm::mat4 m_view = glm::mat4(1);
};

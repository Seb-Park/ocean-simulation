#pragma once

#include "Game/Application.h"
#include "Graphics/global.h"

#include <GLFW/glfw3.h>

class Core
{
public:
    Core();
    ~Core();
    void update(double deltaTime);
    void draw();
    void keyEvent(int key, int action);
    void mousePosEvent(double xpos, double ypos);
    void mouseButtonEvent(int button, int action);
    void scrollEvent(double distance);
    void windowResizeEvent(int width, int height);
    void framebufferResizeEvent(int width, int height);

private:
    std::shared_ptr<Camera> camera;
    std::shared_ptr<Shape> shape;
    std::shared_ptr<ModelTransform> modelTransform;

    bool w_down, a_down, s_down, d_down;
    bool right_mouse_down = false;
    glm::vec2 prev_mouse_pos;
    float m_velocity = .05f;
    Application m_app;
};

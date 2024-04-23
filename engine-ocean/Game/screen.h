#ifndef SCREEN_H
#define SCREEN_H


class Screen
{
public:
    Screen();

    virtual void update(double deltaTime) = 0;
    virtual void draw() = 0;
    virtual void keyEvent(int key, int action) = 0;
    virtual void mousePosEvent(double xpos, double ypos) = 0;
    virtual void mouseButtonEvent(int button, int action) = 0;
    virtual void scrollEvent(double distance) = 0;
    virtual void windowResizeEvent(int width, int height) = 0;
    virtual void framebufferResizeEvent(int width, int height) = 0;
};

#endif // SCREEN_H

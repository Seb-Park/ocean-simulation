#ifndef APPLICATION_H
#define APPLICATION_H

#include "Game/gameplayscreen.h"
#include "Game/screen.h"
#include "Graphics/global.h"

#include <GLFW/glfw3.h>


class Application
{
public:
    Application();
    void update(double deltaTime);
    void draw();
    void keyEvent(int key, int action);
    void mousePosEvent(double xpos, double ypos);
    void mouseButtonEvent(int button, int action);
    void scrollEvent(double distance);
    void windowResizeEvent(int width, int height);
    void framebufferResizeEvent(int width, int height);

private:

    void createKeyInput(int inputVal);
    void inactivateOtherKeys(int keyVal);

//    Screen activeScreen;
//    Screen *game;
//    Screen *menu;
    std::shared_ptr<Screen> activeScreen;
    std::shared_ptr<Screen> game;
    std::shared_ptr<Screen> menu;

    std::map<int, Input> m_input_map;
};


#endif // APPLICATION_H

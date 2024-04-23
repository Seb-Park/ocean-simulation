#ifndef MENUSCREEN_H
#define MENUSCREEN_H
#include "Game/GameWorld.h"
#include "screen.h"
#include <map>


class MenuScreen: public Screen
{
public:
    MenuScreen(std::map<int, Input>& input_map);
    ~MenuScreen();


    void update(double deltaTime) override;
    void draw() override;
    void keyEvent(int key, int action) override;
    void mousePosEvent(double xpos, double ypos) override;
    void mouseButtonEvent(int button, int action) override;
    void scrollEvent(double distance) override;
    void windowResizeEvent(int width, int height) override;
    void framebufferResizeEvent(int width, int height) override;
private:
    std::map<int, Input>& m_input_map;
};

#endif // MENUSCREEN_H

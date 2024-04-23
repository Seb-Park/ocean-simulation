#ifndef GAMEPLAYSCREEN_H
#define GAMEPLAYSCREEN_H
#include "Game/GameWorld.h"
#include "Graphics/global.h"
#include "Graphics/camera.h"
#include <GLFW/glfw3.h>
#include "screen.h"



class GameplayScreen: public Screen
{
public:
    GameplayScreen(std::map<int, Input>& input_map);


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



    GameWorld game;
};

#endif // GAMEPLAYSCREEN_H

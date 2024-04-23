#ifndef DRAWSYSTEM_H
#define DRAWSYSTEM_H


#include "Game/GameObjects/GameObject.h"
#include <vector>
#include "system.h"

class DrawSystem : public System
{
public:
    DrawSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;


private:
    std::map<std::string, std::shared_ptr<GameObject>>& m_gameobjects;
};

#endif // DRAWSYSTEM_H

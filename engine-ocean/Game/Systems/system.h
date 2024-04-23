#ifndef SYSTEM_H
#define SYSTEM_H


#include "Game/GameObjects/GameObject.h"
#include <map>
class System
{
public:
    System();
    virtual void draw() = 0;
    virtual void update(double deltaTime) = 0;
    virtual void scrollEvent(double distance) = 0;
    virtual void mousePosEvent(double xpos, double ypos) = 0;
};

#endif // SYSTEM_H

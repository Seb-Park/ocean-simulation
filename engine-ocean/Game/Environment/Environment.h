#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

class Environment
{
public:
    Environment();
    virtual void draw() = 0;
    virtual void update(double deltaTime) = 0;
};

#endif // ENVIRONMENT_H

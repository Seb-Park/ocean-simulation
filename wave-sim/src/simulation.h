#pragma once

#include "graphics/shape.h"
#include "system.h"

class Shader;

class Simulation
{
public:
    Simulation();

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    void toggleWire();

    void toggleForceRender();

    Eigen::Vector3d calculateFaceNormal(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c);

    bool calculatePointBehindNormal(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, Eigen::Vector3d p);

    System mainSystem;

    enum EstimationMode { EULER, MIDPOINT, ADAPTIVE };

private:
    Shape m_shape;

    Shape m_ground;

    Shape m_extra;

    EstimationMode estimationMode = MIDPOINT;

    void initGround();
    void initExtra();
};

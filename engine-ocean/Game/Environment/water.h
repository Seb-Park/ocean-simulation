#ifndef WATER_H
#define WATER_H
#include "Graphics/global.h"


#include "Environment.h"

class Water : public Environment
{
public:
    Water(std::shared_ptr<Camera> camera);
    void draw() override;
    void update(double deltaTime) override;



private:

    GLuint du_dv_map;
    void initializeQuad();





    std::shared_ptr<Camera> m_camera;
     std::vector<glm::vec3> m_waterPlane;
     std::shared_ptr<ModelTransform> m_waterMT;

     const float WAVE_SPEED = .008f;
     float m_moveFactor = 0.f;




};

#endif // WATER_H

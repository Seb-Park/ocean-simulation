#ifndef GRASSENVIRONMENT_H
#define GRASSENVIRONMENT_H
#include "Game/GameObjects/GameObject.h"
#include "Graphics/global.h"
#include "Graphics/camera.h"
#include <GLFW/glfw3.h>
#include "Environment.h"

class GrassEnvironment : public Environment
{
public:
    GrassEnvironment(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                     std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                     std::shared_ptr<Camera> camera);

    void draw() override;
    void update(double deltaTime) override;
private:
    GLuint VAO, VBO;
    GLuint grass_tex, wind_tex, meadow_tex;
    std::vector<glm::vec3> positions;



    std::map<std::string, std::shared_ptr<GameObject>>& m_rigid_gameobjects;
    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;
    std::shared_ptr<Camera> m_camera;


    void initializeGrassVBO();
    glm::vec2 getRandomOffset(int ub, int numRows);



};

#endif // GRASSENVIRONMENT_H

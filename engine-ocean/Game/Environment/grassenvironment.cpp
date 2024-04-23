#include "grassenvironment.h"
#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/TransformComponent.h"
#include "Game/GameObjects/GameObject.h"
#include "Graphics/global.h"

GrassEnvironment::GrassEnvironment(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                                   std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                   std::shared_ptr<Camera> camera):
    m_rigid_gameobjects(rigid_gameobjects),
    m_dynamic_gameobjects(dynamic_gameobjects),
    m_camera(camera)
{

    initializeGrassVBO();
}

void GrassEnvironment::initializeGrassVBO(){
//    grass_tex = Global::graphics.loadTextureFromFile("/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/meadow_texture_atlas.png").textureID;
//    wind_tex = Global::graphics.loadTextureFromFile("/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/wind_flowmap1.png").textureID;
//    meadow_tex = Global::graphics.loadTextureFromFile("/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/meadow_flower_map.png").textureID;


//    std::vector<glm::vec3> surface_points = m_rigid_gameobjects["ground"]->getComponent<CollisionComponent>()->getCollisionShape<BoundingTriangle>()->getSurfacePoints();
//    positions.insert(std::end(positions), std::begin(surface_points), std::end(surface_points));


//     glGenVertexArrays(1, &VAO);
//     glGenBuffers(1, &VBO);
//     glBindVertexArray(VAO);

//     glBindBuffer(GL_ARRAY_BUFFER, VBO);
//     glBufferData(GL_ARRAY_BUFFER, positions.size()*sizeof(glm::vec3), positions.data(), GL_STATIC_DRAW);

//     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);

//     // unbind
//     glEnableVertexAttribArray(0);
//     glBindBuffer(GL_ARRAY_BUFFER, 0);
}

glm::vec2 GrassEnvironment::getRandomOffset(int ub, int numRows){
    int index = floor(Global::graphics.generateRandomNumbers(0, ub));
    int column = index % numRows;
    float xoffset = static_cast<float>(column)/static_cast<float>(numRows);

    int row = floor(index/numRows);
    float yOffset = static_cast<float>(row)/static_cast<float>(numRows);
    return glm::vec2(xoffset, yOffset);
}

void GrassEnvironment::draw(){
    // draw grass
    glPointSize(8.f);
    Global::graphics.bindShader("grass");

    // activate texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, grass_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, wind_tex);

    glActiveTexture(GL_TEXTURE7);
    glBindTexture(GL_TEXTURE_2D, meadow_tex);

    // update time
    glUniform1f(glGetUniformLocation(Global::graphics.getShaderID("grass"), "u_time"), glfwGetTime()*.1f);
    // update player pos
    glm::vec3 playerPos = m_dynamic_gameobjects.at("player")->getComponent<TransformComponent>()->getPos();
    glUniform3f(glGetUniformLocation(Global::graphics.getShaderID("grass"), "playerPos"), playerPos.x, playerPos.y, playerPos.z);

    //glm::vec2 offset = getRandomOffset(4, 2);
    // update offset
    //glUniform2f(glGetUniformLocation(Global::graphics.getShaderID("grass"), "atlas_offset"), offset.x, offset.y);


    Global::graphics.setCameraData(m_camera);
    Global::graphics.setGlobalData(glm::vec3(.5f));
    glBindVertexArray(VAO);
    glDrawArrays(GL_POINTS, 0, positions.size());
};

void GrassEnvironment::update(double deltaTime){};

#include "water.h"
#include "GLFW/glfw3.h"

Water::Water(std::shared_ptr<Camera> camera):
    m_camera(camera)
{

    initializeQuad();

}


void Water::initializeQuad(){
    m_waterPlane = Global::graphics.addShape("water", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/howl_field_water.obj");
    m_waterMT = std::make_shared<ModelTransform>();
    m_waterMT->setScale(glm::vec3(1.f));
    m_waterMT->setPos(glm::vec3(0.f));
    Global::graphics.setWaterHeight(m_waterMT->getPos().y);

    du_dv_map = Global::graphics.loadTextureFromFile_Repeat("/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/waterDUDV.png").textureID;


}

void Water::draw(){
    Global::graphics.bindShader("water");

    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, Global::graphics.getReflectionTexture());
    glActiveTexture(GL_TEXTURE9);
    glBindTexture(GL_TEXTURE_2D, Global::graphics.getRefractionTexture());
    glActiveTexture(GL_TEXTURE10);
    glBindTexture(GL_TEXTURE_2D, du_dv_map);
    glUniform1f(glGetUniformLocation(Global::graphics.getShaderID("water"), "moveFactor"), m_moveFactor);

    Global::graphics.setCameraData(m_camera);

    Global::graphics.drawShape(Global::graphics.getShape("water"), m_waterMT);

}

void Water::update(double deltaTime){
    m_moveFactor += WAVE_SPEED*deltaTime*.01f;
    // loop back to 0 if it gets to one
    if (m_moveFactor == 1) m_moveFactor = 0;
}



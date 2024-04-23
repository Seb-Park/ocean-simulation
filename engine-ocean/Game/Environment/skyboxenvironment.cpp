#include "skyboxenvironment.h"
#include "Graphics/global.h"


SkyboxEnvironment::SkyboxEnvironment(std::shared_ptr<Camera> camera):
    m_camera(camera),
    m_rotation_mat(std::make_shared<ModelTransform>())
{
    initializeVAO();
}

void SkyboxEnvironment::initializeVAO(){
    skybox_tex = Global::graphics.loadCubeMap(m_skyboxTextureFiles);

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, m_vertices.size()*sizeof(float), m_vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);

    // unbind
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void SkyboxEnvironment::draw(){
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_CULL_FACE);

    Global::graphics.bindShader("skybox");
    //Global::graphics.setCameraData(m_camera);

    // activate texture
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_CUBE_MAP, skybox_tex);

    // manually set view and projection, for non-translating view
    glm::mat4 projection = m_camera->getProjection();
    glm::mat4 view = m_camera->getView();
    view[3] = glm::vec4(0.f);

    glUniformMatrix4fv(glGetUniformLocation(Global::graphics.getShaderID("skybox"), "view"), 1, GL_FALSE, glm::value_ptr(view[0]));
    glUniformMatrix4fv(glGetUniformLocation(Global::graphics.getShaderID("skybox"), "projection"), 1, GL_FALSE, glm::value_ptr(projection[0]));

    // apply rotation matrix
    glUniformMatrix4fv(glGetUniformLocation(Global::graphics.getShaderID("skybox"), "rotation"), 1, GL_FALSE, glm::value_ptr(m_rotation_mat->getRotation()[0]));

    //Global::graphics.setGlobalData(glm::vec3(.5f));
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, m_vertices.size());

    glDepthFunc(GL_LESS);
}

void SkyboxEnvironment::update(double deltaTime){
    m_rotation += ROTATE_SPEED * deltaTime;
    m_rotation_mat->rotate(glm::radians(ROTATE_SPEED), glm::vec3(0.f, 1.f, 0.f));


}

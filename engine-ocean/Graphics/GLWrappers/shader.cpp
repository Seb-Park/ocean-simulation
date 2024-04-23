#include "shader.h"

Shader::Shader(std::vector<GLenum> shaderTypes, std::vector<const char*> filepaths){
    m_handle = ShaderLoader::createShaderProgram(shaderTypes, filepaths);
}

Shader::~Shader(){
    glDeleteProgram(m_handle);
}

void Shader::bind(){
    glUseProgram(m_handle);
}

void Shader::unbind(){
    glUseProgram(0);
}

GLuint Shader::getHandle(){
    return m_handle;
}

void Shader::setMaterial(std::shared_ptr<Material> material){
    ColorSource color_source = material->getColorSource();
    switch(color_source){
        case ColorSource::SOLID_COLOR:
            glUniform1i(glGetUniformLocation(m_handle, "colorSource"), 0);
            glUniform3f(glGetUniformLocation(m_handle, "objColor"), material->getColor().r, material->getColor().g, material->getColor().b);
            break;
        case ColorSource::TEXTURE_COLOR:
            glUniform1i(glGetUniformLocation(m_handle, "colorSource"), 1);
            material->getTexture()->bind();
            glUniform1i(glGetUniformLocation(m_handle, "objTexture"), material->getTexture()->getTexUnitUint());
            break;
        case ColorSource::PER_VERTEX_COLOR:
            glUniform1i(glGetUniformLocation(m_handle, "colorSource"), 2);
            break;
    }
    float shininess = material->getShininess();
    glUniform1f(glGetUniformLocation(m_handle, "shininess"), shininess);
}

void Shader::setCamera(std::shared_ptr<Camera> camera){
    glUniformMatrix4fv(glGetUniformLocation(m_handle, "view"), 1, GL_FALSE, glm::value_ptr(camera->getView()[0]));
    glUniformMatrix4fv(glGetUniformLocation(m_handle, "projection"), 1, GL_FALSE, glm::value_ptr(camera->getProjection()[0]));
    glUniform3f(glGetUniformLocation(m_handle, "worldSpace_camPos"), camera->getPos().x, camera->getPos().y, camera->getPos().z);
    glUniform3f(glGetUniformLocation(m_handle, "skyColor"), 1.f, 1.f, 1.f);

}

void Shader::setModelTransform(std::shared_ptr<ModelTransform> modelTransform){
    glUniformMatrix4fv(glGetUniformLocation(m_handle, "model"), 1, GL_FALSE, glm::value_ptr(modelTransform->getModelMatrix()[0]));
}

void Shader::setModelTransform(glm::mat4 modelMatrix){
    glUniformMatrix4fv(glGetUniformLocation(m_handle, "model"), 1, GL_FALSE, glm::value_ptr(modelMatrix));
}

void Shader::setGlobalCoeffs(glm::vec3 coeffs){
    glUniform3f(glGetUniformLocation(m_handle, "coeffs"), coeffs.x, coeffs.y, coeffs.z);
}

void Shader::setLights(std::vector<std::shared_ptr<Light>> lights){
    int numLights = std::min(int(lights.size()), 16);
    std::vector<int> lightType;
    std::vector<float> lightColor;
    std::vector<float> lightFunction;
    std::vector<float> worldSpace_lightPos;
    std::vector<float> worldSpace_lightDir;
    lightType.resize(numLights);
    lightColor.resize(numLights*3);
    lightFunction.resize(numLights*3);
    worldSpace_lightPos.resize(numLights*3);
    worldSpace_lightDir.resize(numLights*3);
    for(int i = 0; i<numLights; i++){
        lightColor[i*3] = lights[i]->getColor().r;
        lightColor[i*3+1] = lights[i]->getColor().g;
        lightColor[i*3+2] = lights[i]->getColor().b;
        glm::vec3 camLightData;
        switch(lights[i]->getType()){
            case LightType::POINT:
                lightType[i] = 0;
                lightFunction[i*3] = lights[i]->getAttenuation().x;
                lightFunction[i*3+1] = lights[i]->getAttenuation().y;
                lightFunction[i*3+2] = lights[i]->getAttenuation().z;
                worldSpace_lightPos[i*3] = lights[i]->getPos().x;
                worldSpace_lightPos[i*3+1] = lights[i]->getPos().y;
                worldSpace_lightPos[i*3+2] = lights[i]->getPos().z;
                break;
            case LightType::DIRECTIONAL:
                lightType[i] = 1;
                worldSpace_lightDir[i*3] = lights[i]->getDir().x;
                worldSpace_lightDir[i*3+1] = lights[i]->getDir().y;
                worldSpace_lightDir[i*3+2] = lights[i]->getDir().z;
                break;
        }
    }
    glUniform1i(glGetUniformLocation(m_handle, "numLights"), numLights);
    Debug::checkGLError();
    glUniform1iv(glGetUniformLocation(m_handle, "lightType"), numLights, lightType.data());
    Debug::checkGLError();
    glUniform3fv(glGetUniformLocation(m_handle, "lightColor"), numLights, lightColor.data());
    Debug::checkGLError();
    glUniform3fv(glGetUniformLocation(m_handle, "lightFunction"), numLights, lightFunction.data());
    Debug::checkGLError();
    glUniform3fv(glGetUniformLocation(m_handle, "worldSpace_lightPos"), numLights, worldSpace_lightPos.data());
    Debug::checkGLError();
    std::cout<<worldSpace_lightDir[0]<<std::endl;
    glUniform3fv(glGetUniformLocation(m_handle, "worldSpace_lightDir"), numLights, worldSpace_lightDir.data());
    Debug::checkGLError();
}

void Shader::clearLights(){
    glUniform1i(glGetUniformLocation(m_handle, "numLights"), 0);
}

void Shader::setTextUniforms(float screenWidth, float screenHeight, glm::vec3 color){
    glm::mat4 projection = glm::ortho(0.0f, screenWidth, 0.0f, screenHeight);
    glUniformMatrix4fv(glGetUniformLocation(m_handle, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    glUniform3f(glGetUniformLocation(m_handle, "textColor"), color.r, color.g, color.b);
}

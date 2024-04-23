#include "environmentsystem.h"
#include "Game/Environment/grassenvironment.h"
#include "Game/Environment/skyboxenvironment.h"
#include "Game/Environment/water.h"

EnvironmentSystem::EnvironmentSystem(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                                     std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                     std::shared_ptr<Camera> camera)
{
    addEnviron<GrassEnvironment>(std::make_shared<GrassEnvironment>(rigid_gameobjects, dynamic_gameobjects, camera));
    addEnviron<SkyboxEnvironment>(std::make_shared<SkyboxEnvironment>(camera));
    addEnviron<Water>(std::make_shared<Water>(camera));



}

void EnvironmentSystem::draw(){
    for (auto &environ : m_environs){
        environ.second->draw();
    }

}

void EnvironmentSystem::drawNonWater(){
    getEnviron<GrassEnvironment>()->draw();
    getEnviron<SkyboxEnvironment>()->draw();


}

void EnvironmentSystem::update(double deltaTime){
    for (auto &environ : m_environs){
        environ.second->update(deltaTime);
    }

}

void EnvironmentSystem::scrollEvent(double distance){}
void EnvironmentSystem::mousePosEvent(double xpos, double ypos) {}

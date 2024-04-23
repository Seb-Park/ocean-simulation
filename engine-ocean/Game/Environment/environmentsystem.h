#ifndef ENVIRONMENTSYSTEM_H
#define ENVIRONMENTSYSTEM_H
#include "Game/Environment/Environment.h"
#include "Game/Systems/system.h"

class EnvironmentSystem : public System
{
public:
    EnvironmentSystem(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                      std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                      std::shared_ptr<Camera> camera);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;
    void drawNonWater();


    template <typename T>
    void addEnviron(std::shared_ptr<T> &&environ){
        m_environs.put<T>(std::forward<std::shared_ptr<T>>(environ));
    }


    template <class T>
    T* getEnviron(){
        auto comp = m_environs.find<T>();
        assert(comp != m_environs.end());
        return static_cast<T*>(comp->second.get());
    }

    template <class T>
    void removeEnviron(){
        m_environs.remove<T>();
    }
private:
    TypeMap<std::shared_ptr<Environment>> m_environs;

};

#endif // ENVIRONMENTSYSTEM_H

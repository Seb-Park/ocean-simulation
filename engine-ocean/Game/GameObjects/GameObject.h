#ifndef GAMEOBJECT_H
#define GAMEOBJECT_H


#include "Game/Components/Component.h"
#include "Game/TypeMap.h"
#include <memory>
#include "Graphics/global.h"

class GameObject
{
public:
    GameObject();
    ~GameObject(){
        std::cout << "deleting" << std::endl;
    }

    template <typename T>
    void addComponent(std::unique_ptr<T> &&component){
        m_components.put<T>(std::forward<std::unique_ptr<T>>(component));
    }

    template <typename T>
    bool hasComponent(){
        return m_components.contains<T>();
    }

    template <class T>
    T* getComponent(){
        auto comp = m_components.find<T>();
        assert(comp != m_components.end());
        return static_cast<T*>(comp->second.get());
    }

    template <class T>
    void removeComponent(){
//        auto comp = m_components.find<T>();
//        assert(comp != m_components.end());
        m_components.remove<T>();
    }



private:
    void setPlayerPos();
    void makePlayer();

    TypeMap<std::unique_ptr<Component>> m_components;
    std::shared_ptr<Shape> m_object;
    std::shared_ptr<ModelTransform> m_object_mt;

};

#endif // GAMEOBJECT_H

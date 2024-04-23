#ifndef UIELEMENT_H
#define UIELEMENT_H


#include "Game/Systems/UI/UITextures/uitexture.h"
#include "Game/TypeMap.h"
#include <memory>
class UIElement
{
public:
    UIElement();

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
        m_components.remove<T>();
    }

private:

     TypeMap<std::unique_ptr<UITexture>> m_components;
};

#endif // UIELEMENT_H

#ifndef INVENTORYITEM_H
#define INVENTORYITEM_H


#include "Game/Systems/UI/UITextures/uidisplay.h"
class InventoryItem
{
public:
    InventoryItem(const char* worldLabelFile, const char* inventoryLabelFile,
                  glm::vec2 worldLabelScale,
                  glm::vec2 inventoryLabelScale,
                  std::set<std::string>& m_shownScreens);
    void updateWorldLabelPos(glm::vec2 pos);
    glm::mat4 getWorldLabelTransformationMat();
    GLuint getWorldLabelTexID();
    float getWorldLabelTexAspect();
    std::shared_ptr<UIDisplay> getUIDisplay();





private:
    std::shared_ptr<UIDisplay> m_world_label;
    std::shared_ptr<UIDisplay> m_inventory_label;

    int m_inventoryCount = 0;

};

#endif // INVENTORYITEM_H

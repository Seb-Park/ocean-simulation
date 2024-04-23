#include "inventoryitem.h"

InventoryItem::InventoryItem(const char* worldLabelFile, const char* inventoryLabelFile,
                             glm::vec2 worldLabelScale,
                             glm::vec2 inventoryLabelScale,
                             std::set<std::string>& m_shownScreens)
{
    TextureData worldLabelTex = Global::graphics.loadTextureFromFile(worldLabelFile);
    TextureData inventoryLabelTex = Global::graphics.loadTextureFromFile(inventoryLabelFile);

    m_world_label = std::make_shared<UIDisplay>(worldLabelTex, glm::vec2(0.f), worldLabelScale, m_shownScreens, AspectRatio::LAND_FIT);
    m_inventory_label = std::make_shared<UIDisplay>(inventoryLabelTex, glm::vec2(0.f), inventoryLabelScale, m_shownScreens, AspectRatio::LAND_FIT);
}

void InventoryItem::updateWorldLabelPos(glm::vec2 pos){
    m_world_label->setPos(pos);
}

glm::mat4 InventoryItem::getWorldLabelTransformationMat(){
    return m_world_label->getTransformationMat();
}

GLuint InventoryItem::getWorldLabelTexID(){
    return m_world_label->getTexID();
}

float InventoryItem::getWorldLabelTexAspect(){
    return m_world_label->getTextureScaleAspect();
}

std::shared_ptr<UIDisplay> InventoryItem::getUIDisplay(){
    return m_world_label;
}

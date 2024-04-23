#ifndef INVENTORYSYSTEM_H
#define INVENTORYSYSTEM_H


#include "Game/GameWorld.h"
#include <map>
#include "inventoryitem.h"

//struct InventoryItem{
//    TextureData lootLabelUI;
//    TextureData lootInventoryUI;
//    int inventoryCount = 0;

//    InventoryItem(TextureData label, TextureData inventory):
//        lootLabelUI(label),
//        lootInventoryUI(inventory)
//    {}
//};

class InventorySystem : public System
{
public:
    InventorySystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                    std::map<std::string, std::shared_ptr<GameObject>>& all_gameobjects,
                    std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables,
                                     std::map<int, Input>& input_map,
                                        std::shared_ptr<Camera>& camera,
                    std::set<std::string>& m_shownScreens);
    ~InventorySystem();
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;
    void drawUIText();
    void onWindowResize(int width, int height);



private:
    void addToInventory(std::shared_ptr<GameObject> &go);
    bool withinPlayerReach(const std::shared_ptr<GameObject> &go);
    void initializeInventory(std::string itemName, const char* label_filename, const char* inventory_filename);





    int m_bamboo_count = 0;
    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;

    // maps loot id to game objects/associated textures
    std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& m_lootables; // holds position
    std::map<std::string, std::shared_ptr<InventoryItem>> m_inventoryItems;

    std::map<int, Input>& m_input_map;
     std::shared_ptr<Camera>& m_camera;

     //std::set<std::string> m_inventoryScreens; // stores ids of all scenes in which inventory is rendered
     std::set<std::string>& m_shownScreens;

     std::vector<float> m_quadPos = {
         -1.0f, 1.0f,
         -1.0f, -1.0f,
         1.0f, 1.0f,
         1.0f, -1.0f
     };

     GLuint m_screenVAO;
     GLuint m_texID;

     glm::vec2 m_offset = glm::vec2(.2f);

     std::shared_ptr<UIDisplay> m_sparkle;

     std::map<std::string, std::shared_ptr<GameObject>>& m_all_gameobjects;
};

#endif // INVENTORYSYSTEM_H

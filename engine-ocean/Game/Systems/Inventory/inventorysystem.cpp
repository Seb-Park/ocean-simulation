#include "inventorysystem.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/drawcomponent.h"
#include "Game/Components/transformcomponent.h"
#include "Game/GameObjects/GameObject.h"
#include "Game/GameWorld.h"
#include <memory>

InventorySystem::InventorySystem(std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                 std::map<std::string, std::shared_ptr<GameObject>>& all_gameobjects,
                                 std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables,
                                 std::map<int, Input>& input_map,
                                  std::shared_ptr<Camera>& camera,
                                 std::set<std::string>& m_shownScreens):
    m_dynamic_gameobjects(dynamic_gameobjects),
    m_input_map(input_map),
    m_camera(camera),
    m_lootables(lootables),
    m_shownScreens(m_shownScreens),
    m_all_gameobjects(all_gameobjects)
{
    int m_bamboo_count = 0;
    m_screenVAO = Global::graphics.makeVAO(m_quadPos);

    //m_texID = Global::graphics.loadTextureFromFile("/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/enterbutton.png").textureID;


    initializeInventory("mushroom",
                        "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/mushroom_loot.png",
                        "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/mushroom_loot.png");
    initializeInventory("amantia",
                        "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/amantia_loot.png",
                        "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/amantia_loot.png");

    TextureData sparkle = Global::graphics.loadTextureFromFile("/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/sparkle.png");
    m_sparkle = std::make_shared<UIDisplay>(sparkle, glm::vec2(0.f), glm::vec2(.05f), m_shownScreens, AspectRatio::LAND_FIT);

}

InventorySystem::~InventorySystem(){
    glDeleteVertexArrays(1, &m_screenVAO);
}

void InventorySystem::initializeInventory(std::string itemName, const char* label_filename, const char* inventory_filename){
    std::shared_ptr<InventoryItem> inventoryItem = std::make_shared<InventoryItem>(label_filename, inventory_filename, glm::vec2(.2f), glm::vec2(.5f), m_shownScreens);
    m_inventoryItems.insert({itemName, inventoryItem});
}


void InventorySystem::addToInventory(std::shared_ptr<GameObject> &go){
    // if G is pressed down, remove go from screen and add it to the inventory
    if (m_input_map.at(GLFW_KEY_P).isActive){
        if (go->hasComponent<DrawComponent>()){
            go->removeComponent<DrawComponent>();
            go->removeComponent<CollisionComponent>();
            m_bamboo_count += 1;
        }
    }
}

bool InventorySystem::withinPlayerReach(const std::shared_ptr<GameObject> &go){
        glm::vec3 player_pos = m_dynamic_gameobjects.at("player")->getComponent<TransformComponent>()->getPos();

        // be careful --> the below line will crash if its attempting to get an object with multiple
        // transforms, like the ground!
        glm::vec3 go_pos = go->getComponent<TransformComponent>()->getPos();

        float vicinity = 3.f;
        // is player and game object are close in either x or z dir
        if ((abs(player_pos.x - go_pos.x) < vicinity)
             && (abs(player_pos.z - go_pos.z) < vicinity)
             && (abs(player_pos.y - go_pos.y) < vicinity)){
            return true;
        }
    return false;
}

void InventorySystem::drawUIText(){
    Global::graphics.drawUIText(Global::graphics.getFont("opensans"), std::to_string(m_bamboo_count) + " bamboo collected.", glm::ivec2(20, 90), AnchorPoint::TopLeft, Global::graphics.getFramebufferSize().x, 0.2f, 0.1f, glm::vec3(1, 1, 1));
    Global::graphics.drawUIText(Global::graphics.getFont("opensans"), "Collect 20 bamboos to clear the game.", glm::ivec2(20, 110), AnchorPoint::TopLeft, Global::graphics.getFramebufferSize().x, 0.2f, 0.1f, glm::vec3(1, 1, 1));

    if (m_bamboo_count  >= 20){
        Global::graphics.drawUIText(Global::graphics.getFont("opensans"), "You Won! Press 'R' to restart.", glm::ivec2(200, 300), AnchorPoint::TopLeft, Global::graphics.getFramebufferSize().x, 0.4f, 0.4f, glm::vec3(1, 1, 1));
    }
}




void InventorySystem::draw(){
    // bind shader
    Global::graphics.bindShader("inventory");

    // for window resizing
    Global::graphics.setCameraData(m_camera);

    // enable alpha blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // disable depth testing
    glDisable(GL_DEPTH_TEST);

    // for every item in inventory, compute screen coordinates and draw if near player and on the screen
    glm::mat4 proj = m_camera->getProjection();
    glm::mat4 view = m_camera->getView();

    // if home page is rendered
    if (m_shownScreens.contains("home")){
        glBindVertexArray(m_screenVAO);
        glEnableVertexAttribArray(0);

        for (auto &lootID : m_lootables){
            for (auto &lootObject : lootID.second){
                glm::vec3 pos = lootObject->getComponent<TransformComponent>()->getMT()->getPos();
                glm::vec4 projected = proj*view*glm::vec4(pos, 1.f);
                glm::vec2 screen = glm::vec2(projected.x/projected.w, projected.y/projected.w);

                // if visible on screen
                if ((-1.f < screen.x && screen.x < 1.f) && (-1.f < screen.y && screen.y < 1.f)){



                    // update texture pos, and draw associated texture
                    if (withinPlayerReach(lootObject)){
                        // update pos
                        m_inventoryItems[lootID.first]->updateWorldLabelPos(screen + m_offset);

                        // get uniforms
                        glm::mat4 transMat = m_inventoryItems[lootID.first]->getWorldLabelTransformationMat();
                        GLuint texID = m_inventoryItems[lootID.first]->getWorldLabelTexID();
                        float texAspect = m_inventoryItems[lootID.first]->getWorldLabelTexAspect();

                        glActiveTexture(GL_TEXTURE6);
                        glBindTexture(GL_TEXTURE_2D, texID);
                        glm::vec2 texScale = glm::vec2(1.f, texAspect);
                        glUniform2f(glGetUniformLocation(Global::graphics.getShaderID("inventory"), "textureScale"), texScale.x, texScale.y);
                        glUniformMatrix4fv(glGetUniformLocation(Global::graphics.getShaderID("inventory"), "transform"), 1, GL_FALSE, glm::value_ptr(transMat[0]));
                        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

                    }
                    // draw sparkle
                    else {
                        // update pos
                        m_sparkle->setPos(screen);

                        // get uniforms
                        glm::mat4 transMat = m_sparkle->getTransformationMat();
                        GLuint texID = m_sparkle->getTexID();
                        float texAspect = m_sparkle->getTextureScaleAspect();

                        glActiveTexture(GL_TEXTURE6);
                        glBindTexture(GL_TEXTURE_2D, texID);
                        glm::vec2 texScale = glm::vec2(1.f/texAspect, 1.f);
                        glUniform2f(glGetUniformLocation(Global::graphics.getShaderID("inventory"), "textureScale"), texScale.x, texScale.y);
                        glUniformMatrix4fv(glGetUniformLocation(Global::graphics.getShaderID("inventory"), "transform"), 1, GL_FALSE, glm::value_ptr(transMat[0]));
                        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

                    }
                }

            }

        }

        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
    }
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}


void InventorySystem::update(double deltaTime){

    // remove from m_lootables to remove labels and drawing
    for (auto &lootGroup : m_lootables){
        for (auto itr = lootGroup.second.begin(); itr != lootGroup.second.end();){
            auto& loot = *itr;
            if (withinPlayerReach(loot) && m_input_map.at(GLFW_KEY_P).isActive){
                loot->removeComponent<DrawComponent>();
                itr = lootGroup.second.erase(itr);
            } else {
                ++itr;
            }
         }
    }

    // remove game object entirely based on if its not drawn
//    for(auto it = m_all_gameobjects.begin(); it != m_all_gameobjects.end();){
//        if (!it->second->hasComponent<DrawComponent>()){
//            it = m_all_gameobjects.erase(it);
//        }
//    }
}

void InventorySystem::onWindowResize(int width, int height){
    for (auto &lootID : m_lootables){
        m_inventoryItems[lootID.first]->getUIDisplay()->setWindowPos(width, height);

    }
}

void InventorySystem::scrollEvent(double distance){}
void InventorySystem::mousePosEvent(double xpos, double ypos){}


#include "uisystem.h"
#include "Game/GameWorld.h"
#include "Game/Systems/UI/ButtonAction/showwindowaction.h"
#include "Game/Systems/UI/UITextures/UIButton.h"
#include "Game/Systems/UI/UITextures/uidisplay.h"
#define GLFW_POINTING_HAND_CURSOR

UISystem::UISystem(std::shared_ptr<Camera> camera,
                   std::map<int, Input>& input_map,
                   std::set<std::string>& shownScreens):
    m_camera(camera),
    m_input_map(input_map),
    m_shownScreens(shownScreens)
{


    initializeStartScreen();
    initializeScreenMap();
    //m_pointerCursor = glfwCreateStandardCursor(GLFW_POINTING_HAND_CURSOR);

}

UISystem::~UISystem(){
    for (auto &screenID : m_all_screens){
        glDeleteVertexArrays(1, &screenID.second->screenVAOID);
    }
}

void UISystem::initializeScreenMap(){
    initializeProfileScreen();
    initializeInventory();
    initializeScreen();

    m_shownScreens.insert("home");

}

void UISystem::initializeScreen(){
    std::shared_ptr<UIScreen> home = std::make_shared<UIScreen>();
    home->isActive = true;
    home->screenVAOID = Global::graphics.makeVAO(m_quadPos);
    m_quad_numVertices = 4;

    float topLevel = .85f;

    makeButtonElement(home, "profileButton", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/profileicon.png", glm::vec2(.85f,topLevel), glm::vec2(.05f), AspectRatio::LAND_FILL);
    std::shared_ptr<ButtonAction> showProfile = std::make_shared<ShowWindowAction>(m_all_screens, m_shownScreens, "profile");
    home->screenElements.at("profileButton")->getComponent<UIButton>()->addButtonAction(showProfile);

    makeButtonElement(home, "inventoryButton", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/inventoryicon.png", glm::vec2(.75f, topLevel), glm::vec2(.05f), AspectRatio::LAND_FILL);
    std::shared_ptr<ButtonAction> showInv = std::make_shared<ShowWindowAction>(m_all_screens, m_shownScreens, "inventory");
    home->screenElements.at("inventoryButton")->getComponent<UIButton>()->addButtonAction(showInv);

    makeButtonElement(home, "questsButton", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/questicon.png", glm::vec2(.65f, topLevel), glm::vec2(.05f), AspectRatio::LAND_FILL);
    makeButtonElement(home, "settingsButton", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/settings_icon.png", glm::vec2(-.85f, topLevel), glm::vec2(.05f), AspectRatio::LAND_FILL);
    makeDisplayElement(home, "healthbar", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/healthbar.png", glm::vec2(0.f, -.5f), glm::vec2(.4f), AspectRatio::LAND_FIT);

//    TextureData reflectTex;
//    reflectTex.textureID = Global::graphics.getReflectionTexture();
//    reflectTex.width = Global::graphics.REFLECTION_W;
//    reflectTex.height = Global::graphics.REFLECTION_H;
//    makeDisplayElement(home, "reflection", reflectTex, glm::vec2(.5f, -.5f), glm::vec2(.3f), AspectRatio::LAND_FIT);

//    TextureData refractTex;
//    refractTex.textureID = Global::graphics.getRefractionTexture();
//    refractTex.width = Global::graphics.REFRACTION_W;
//    refractTex.height = Global::graphics.REFRACTION_H;
//    makeDisplayElement(home, "refract", refractTex, glm::vec2(-.5f, -.5f), glm::vec2(.3f), AspectRatio::LAND_FIT);




    home->elementsDepthOrder.push_back("profileButton");
    home->elementsDepthOrder.push_back("inventoryButton");
    home->elementsDepthOrder.push_back("questsButton");
    home->elementsDepthOrder.push_back("settingsButton");
    home->elementsDepthOrder.push_back("healthbar");
   // home->elementsDepthOrder.push_back("reflection");
   // home->elementsDepthOrder.push_back("refract");



    m_all_screens.insert({"home", home});
}

void UISystem::initializeProfileScreen(){

    std::shared_ptr<UIScreen> profile = std::make_shared<UIScreen>();
    profile->isActive = false;
    profile->screenVAOID = Global::graphics.makeVAO(m_quadPos);
    m_quad_numVertices = 4;

    float topLevel = .85f;

    // eventually should be display, not button
    makeDisplayElement(profile, "bg", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/inventory_bg.png", glm::vec2(0.f), glm::vec2(1.0f));
    makeDisplayElement(profile, "profileDisplay", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/mouse_profile.png", glm::vec2(0.f), glm::vec2(1.0f), AspectRatio::LAND_FIT);
    makeButtonElement(profile, "closeButton", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/close_icon.png", glm::vec2(0.f), glm::vec2(.05f),  AspectRatio::LAND_FILL, true, CornerPosition::TOPRIGHT);
    profile->screenElements.at("closeButton")->getComponent<UIButton>()->setWindowToClose("profile");

    profile->elementsDepthOrder.push_back("bg");
    profile->elementsDepthOrder.push_back("profileDisplay");
    profile->elementsDepthOrder.push_back("closeButton");

    m_all_screens.insert({"profile", profile});
}

void UISystem::initializeInventory(){
    std::shared_ptr<UIScreen> inv = std::make_shared<UIScreen>();
    inv->isActive = false;
    inv->screenVAOID = Global::graphics.makeVAO(m_quadPos);
    m_quad_numVertices = 4;

    makeDisplayElement(inv, "bg", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/inventory_bg.png", glm::vec2(0.f), glm::vec2(1.0f));
    makeDisplayElement(inv, "inventoryDisplay", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/inventory_page.png", glm::vec2(0.f), glm::vec2(1.0f), AspectRatio::LAND_FIT);
    makeButtonElement(inv, "closeButton", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/close_icon.png", glm::vec2(0.f), glm::vec2(.05f),  AspectRatio::LAND_FILL, true, CornerPosition::TOPRIGHT);
    inv->screenElements.at("closeButton")->getComponent<UIButton>()->setWindowToClose("inventory");
    inv->elementsDepthOrder.push_back("bg");
    inv->elementsDepthOrder.push_back("inventoryDisplay");
    inv->elementsDepthOrder.push_back("closeButton");

    m_all_screens.insert({"inventory", inv});
}

void UISystem::initializeStartScreen(){

    std::shared_ptr<UIScreen> start = std::make_shared<UIScreen>();
    start->isActive = true;
    start->screenVAOID = Global::graphics.makeVAO(m_quadPos);

    float topLevel = .85f;

    // eventually should be display, not button
    makeDisplayElement(start, "startDisplay", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/title.png", glm::vec2(0.f), glm::vec2(1.0f), AspectRatio::LAND_FIT);
    makeButtonElement(start, "startButton", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/enterbutton.png", glm::vec2(0.f, -.57f), glm::vec2(.19f), AspectRatio::LAND_FIT, true);
    start->screenElements.at("startButton")->getComponent<UIButton>()->setWindowToClose("start");

    start->elementsDepthOrder.push_back("startDisplay");
    start->elementsDepthOrder.push_back("startButton");

    m_all_screens.insert({"start", start});
    m_shownScreens.insert("start");
}

glm::vec2 UISystem::drawAspect(AspectRatio aspectType, float textureAspect){
    switch(aspectType){
    case FIT_SCREEN:
        return glm::vec2(1.f);
        break;
    case LAND_FIT:case PORTRAIT_FILL:
        return glm::vec2(1.f, textureAspect);
        break;
    case LAND_FILL: case PORTRAIT_FIT:
        return glm::vec2(1.f/textureAspect, 1.f);
        break;
    }
}



void UISystem::makeButtonElement(std::shared_ptr<UIScreen> &screen,
                                 std::string elementName, const char* filename,
                                 const glm::vec2 pos, const glm::vec2 scale,
                                 AspectRatio aspectRatio,
                                 bool isCloseButton,
                                 CornerPosition corner){
    TextureData tex = Global::graphics.loadTextureFromFile(filename);

    std::shared_ptr<UIElement> uiElement = std::make_shared<UIElement>();

    uiElement->addComponent<UIButton>(std::make_unique<UIButton>(tex, pos, scale, m_shownScreens, isCloseButton, corner, aspectRatio));
    screen->screenElements.insert({elementName, uiElement});
}

void UISystem::makeDisplayElement(std::shared_ptr<UIScreen> &screen,
                                  std::string elementName, const char* filename,
                                  const glm::vec2 pos, const glm::vec2 scale,
                                  AspectRatio aspectRatio){
    TextureData tex = Global::graphics.loadTextureFromFile(filename);

    std::shared_ptr<UIElement> uiElement = std::make_shared<UIElement>();
    uiElement->addComponent<UIDisplay>(std::make_unique<UIDisplay>(tex, pos, scale, m_shownScreens, aspectRatio));
    screen->screenElements.insert({elementName, uiElement});
}

void UISystem::makeDisplayElement(std::shared_ptr<UIScreen> &screen,
                                  std::string elementName, TextureData &tex,
                                  const glm::vec2 pos, const glm::vec2 scale,
                                  AspectRatio aspectRatio){

    std::shared_ptr<UIElement> uiElement = std::make_shared<UIElement>();
    uiElement->addComponent<UIDisplay>(std::make_unique<UIDisplay>(tex, pos, scale, m_shownScreens, aspectRatio));
    screen->screenElements.insert({elementName, uiElement});
}

void UISystem::renderScreen(){
    // bind shader
    Global::graphics.bindShader("ui");

    // for window resizing
    Global::graphics.setCameraData(m_camera);

    // enable alpha blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // disable depth testing
    glDisable(GL_DEPTH_TEST);

    for (auto &screenID : m_shownScreens){
            // bind quad
            glBindVertexArray(m_all_screens[screenID]->screenVAOID);
            glEnableVertexAttribArray(0);

            // loop through inventory ui elements

            // loop through ui elements
            for (auto &gui: m_all_screens[screenID]->elementsDepthOrder){
                if (m_all_screens[screenID]->screenElements[gui]->hasComponent<UIDisplay>()){
                    GLuint texID;
                    if (gui == "reflection"){
                        //std::cout << "rendering reflection" << std::endl;
                        texID = Global::graphics.getReflectionTexture();
                    } else if (gui == "refract"){
                        texID = Global::graphics.getRefractionTexture();
                    } else {
                        texID = m_all_screens[screenID]->screenElements[gui]->getComponent<UIDisplay>()->getTexID();
                    }
                    glm::mat4 transMat = m_all_screens[screenID]->screenElements[gui]->getComponent<UIDisplay>()->getTransformationMat();
                    float texAspect = m_all_screens[screenID]->screenElements[gui]->getComponent<UIDisplay>()->getTextureScaleAspect();
                    AspectRatio ratio = m_all_screens[screenID]->screenElements[gui]->getComponent<UIDisplay>()->getAspectRatio();

                    glActiveTexture(GL_TEXTURE5);
                    glBindTexture(GL_TEXTURE_2D, texID);

                    glm::vec2 texScale = drawAspect(ratio, texAspect);

                    glUniform2f(glGetUniformLocation(Global::graphics.getShaderID("ui"), "textureScale"), texScale.x, texScale.y);
                    glUniformMatrix4fv(glGetUniformLocation(Global::graphics.getShaderID("ui"), "transform"), 1, GL_FALSE, glm::value_ptr(transMat[0]));
                    glDrawArrays(GL_TRIANGLE_STRIP, 0, m_quad_numVertices);
                }
                if (m_all_screens[screenID]->screenElements[gui]->hasComponent<UIButton>()){
                    GLuint texID = m_all_screens[screenID]->screenElements[gui]->getComponent<UIButton>()->getTexID();
                    glm::mat4 transMat = m_all_screens[screenID]->screenElements[gui]->getComponent<UIButton>()->getTransformationMat();
                    float texAspect = m_all_screens[screenID]->screenElements[gui]->getComponent<UIButton>()->getTextureScaleAspect();
                    AspectRatio ratio = m_all_screens[screenID]->screenElements[gui]->getComponent<UIButton>()->getAspectRatio();

                    glActiveTexture(GL_TEXTURE5);
                    glBindTexture(GL_TEXTURE_2D, texID);

                    glm::vec2 texScale = drawAspect(ratio, texAspect);

                    glUniform2f(glGetUniformLocation(Global::graphics.getShaderID("ui"), "textureScale"), texScale.x, texScale.y);
                    glUniformMatrix4fv(glGetUniformLocation(Global::graphics.getShaderID("ui"), "transform"), 1, GL_FALSE, glm::value_ptr(transMat[0]));
                    glDrawArrays(GL_TRIANGLE_STRIP, 0, m_quad_numVertices);


                }

            }
            glDisableVertexAttribArray(0);
            glBindVertexArray(0);
    }
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}



void UISystem::draw(){
    renderScreen();
}

glm::mat4 UISystem::getScaleMatrix(glm::vec2 scale) {
    glm::mat4 M = glm::mat4(1.f);
    M[0][0] = scale.x *  (m_aspectRatio);
    M[1][1] = scale.y;
    M[2][2] = 0.f;
    return M;
}

glm::mat4 UISystem::makeTransformationMat(glm::vec2 translation, glm::vec2 scale){
    glm::mat4 transMat = glm::mat4(1.f);
    transMat[3] = glm::vec4(translation.x, translation.y, 0.f, 1.f);

    glm::mat4 scaleMat = getScaleMatrix(glm::vec2(scale));

    return transMat*scaleMat;
}

void UISystem::update(double deltaTime){
    // allow for one more update cycle before deactivating click
    if (m_input_map[GLFW_MOUSE_BUTTON_LEFT].isClicked){
        // check if clicked inside gui
        // only act on the front-most gui being clicked
        std::set<std::string>::iterator it = m_shownScreens.begin();
        while(it != m_shownScreens.end()) {
            std::set<std::string>::iterator current = it++;
            std::string screenID = *current;
            // render all active screens
            if (m_all_screens[screenID]->isActive){
                // loop through ui elements
                for (auto &gui: m_all_screens[screenID]->screenElements){
                    if (gui.second->hasComponent<UIButton>()){
                        bool isInside = false;
                        Bounds2f buttonBounds = gui.second->getComponent<UIButton>()->getBounds();
                        isInside = checkInsideGUI(m_mousepos, buttonBounds);
                        if (isInside){
                            bool toClose = gui.second->getComponent<UIButton>()->onButtonPress();
                            if (toClose){
                                m_shownScreens.erase(current);

                            }
                        }
                    }
                }
            }
        }

        if (m_input_map[GLFW_MOUSE_BUTTON_LEFT].checkClickTime > 1){
            m_input_map[GLFW_MOUSE_BUTTON_LEFT].isClicked = false;
        } else {
            m_input_map[GLFW_MOUSE_BUTTON_LEFT].checkClickTime ++;
        }
    }
}

bool UISystem::checkInsideGUI(glm::vec2 mouseClickPos, Bounds2f guiBounds){
    if (mouseClickPos.x < guiBounds.max.x
            && mouseClickPos.x > guiBounds.min.x
            && mouseClickPos.y > guiBounds.max.y
            && mouseClickPos.y < guiBounds.min.y){
        return true;
    }

    return false;
}


void UISystem::scrollEvent(double distance){}

void UISystem::mousePosEvent(double xpos, double ypos){
    m_mousepos = glm::vec2(xpos, ypos);

    for (auto &screenID : m_shownScreens){
        Global::graphics.bindShader("ui");
        //glUniform1i(glGetUniformLocation(Global::graphics.getShaderID("ui"), "hovering"), false);
        // render all active screens
        if (m_all_screens[screenID]->isActive){
            // loop through ui elements
            for (auto &gui: m_all_screens[screenID]->screenElements){
                if (gui.second->hasComponent<UIButton>()){
                    bool isInside;
                    Bounds2f buttonBounds = gui.second->getComponent<UIButton>()->getBounds();
                    isInside = checkInsideGUI(m_mousepos, buttonBounds);
                    if (isInside){
                        //glfwSetCursor(window, m_pointerCursor);
                        // hover effect
                        //glUniform1i(glGetUniformLocation(Global::graphics.getShaderID("ui"), "hovering"), true);
                    }
                }
            }
        }
    }
}

void UISystem::onWindowResize(int width, int height){
    m_windowWidth = width;
    m_windowHeight = height;
    m_aspectRatio = static_cast<float>(m_windowHeight)/static_cast<float>(m_windowWidth);

    for (auto &screen : m_all_screens){
            // loop through ui elements
            for (auto &gui: screen.second->screenElements){
                if (gui.second->hasComponent<UIButton>()){
                    gui.second->getComponent<UIButton>()->setWindowPos(width, height);
                }
                if (gui.second->hasComponent<UIDisplay>()){
                    gui.second->getComponent<UIDisplay>()->setWindowPos(width, height);
                }
            }
        }
}

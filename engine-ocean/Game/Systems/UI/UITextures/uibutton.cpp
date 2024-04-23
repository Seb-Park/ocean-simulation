#include "UIButton.h"
#include "Game/Systems/UI/UITextures/UIDisplay.h"
#include <set>

UIButton::UIButton(TextureData tex, glm::vec2 pos, glm::vec2 scale, std::set<std::string>& shownScreens,
                   bool isCloseButton, CornerPosition corner, AspectRatio ratio):
    m_tex(tex),
    m_shownScreens(shownScreens),
    m_aspectRatio(ratio)
{

    // set variables
    m_isCloseButton = isCloseButton;
    m_pos = pos;
    m_scale = scale;
    m_tex_aspectRatio = static_cast<float>(m_tex.height)/static_cast<float>(m_tex.width);

    setToCorner(corner, 640, 480);
    setWindowPos(640, 480);
    setTransformationMat(m_pos, scale);

}
UIButton::~UIButton(){
    glDeleteTextures(1, &m_tex.textureID);
}


AspectRatio UIButton::getAspectRatio(){
    return m_aspectRatio;
}

glm::mat4 UIButton::getScaleMatrix(glm::vec2 scale) {
    glm::mat4 M = glm::mat4(1.f);
    M[0][0] = scale.x; //* (m_screenImageHeight/m_screenImageWidth);
    M[1][1] = scale.y; //* (m_tex_aspectRatio);
    M[2][2] = 1.f;
    return M;
}

void UIButton::setTransformationMat(glm::vec2 translation, glm::vec2 scale){
    glm::mat4 transMat = glm::mat4(1.f);
    transMat[3] = glm::vec4(translation.x, translation.y, 0.f, 1.f);
    glm::mat4 scaleMat = getScaleMatrix(glm::vec2(scale));
    m_transformationMat = transMat*scaleMat;
}

glm::mat4 UIButton::getTransformationMat(){
    return m_transformationMat;
}

float UIButton::getTextureRatio(){
    return m_tex_aspectRatio;
}


void UIButton::draw(){}

GLuint UIButton::getTexID(){
    return m_tex.textureID;
}

glm::vec2 UIButton::getPos(){
 return m_pos;
}

glm::vec2 UIButton::getScale(){
 return m_scale;
}

//void UIButton::setTexID(GLuint &newTexID){
//    m_tex.textureID = newTexID;
//}

void UIButton::setWindowPos(int width, int height){
    m_windowHeight = static_cast<float>(height);
    m_windowWidth = static_cast<float>(width);

    // find where on window it is, for bound checking
    float xpos = .5f*(m_pos.x + 1.f)*m_windowWidth;
    float ypos = (1.f-.5f*(m_pos.y + 1.f))*m_windowHeight;
    m_windowPos = glm::vec2(xpos, ypos);

    // set everything according to window dimensions -- this is for bound checking
    m_toScreenScale = m_windowHeight*m_scale.y;
    m_screenImageHeight = m_toScreenScale;
    m_screenImageWidth = m_toScreenScale * m_tex_aspectRatio;

    float windowRatio = m_windowHeight/m_windowWidth;
    m_textureAspect = windowRatio / m_tex_aspectRatio;

    // calculate window bounds
    glm::vec2 halfDimensions = glm::vec2(m_screenImageWidth, m_screenImageHeight)*.5f;
    m_bounds.max = glm::vec2(m_windowPos.x + halfDimensions.x, m_windowPos.y - halfDimensions.y);
    m_bounds.min = glm::vec2(m_windowPos.x - halfDimensions.x, m_windowPos.y + halfDimensions.y);
}

float UIButton::getTextureScaleAspect(){
    return m_textureAspect;
}

void UIButton::setParentDisplay(std::string parent){
    m_parentDisplay = parent;

}

std::string UIButton::getParentDisplay(){
    return m_parentDisplay;

}

CornerPosition UIButton::getCornerPos(){
    return m_cornerPos;

}

void UIButton::setToCorner(const CornerPosition corner, int width, int height){

    m_toScreenScale = static_cast<float>(height)*m_scale.y;
    m_screenImageHeight = m_toScreenScale;
    m_screenImageWidth = m_toScreenScale * m_tex_aspectRatio;

    // in texture space
    float xtrans = m_screenImageWidth/(.5f*static_cast<float>(width));
    float ytrans = m_screenImageHeight/(.5f*static_cast<float>(height));

    // find where window pos should be
    switch(corner){
        case TOPLEFT:
            m_pos = glm::vec2(-1+xtrans, 1-ytrans);

        break;
        case TOPRIGHT:
            m_pos = glm::vec2(1-xtrans, 1-ytrans);

        break;
        case BOTTOMLEFT:
            m_pos = glm::vec2(-1+xtrans, -1+ytrans);

        break;
        case BOTTOMRIGHT:
            m_pos = glm::vec2(1-xtrans, -1+ytrans);
        break;
        default:
        break;
    }
}

glm::vec2 UIButton::getWindowPos(){
    return m_windowPos;
}

int UIButton::getHeight(){
    return m_tex.height;
}
int UIButton::getWidth(){
    return m_tex.width;
}

int UIButton::getScreenHeight(){
    return m_screenImageHeight;
}
int UIButton::getScreenWidth(){
    return m_screenImageWidth;
}

// remember that origin is top left corner!!
Bounds2f UIButton::getBounds(){

    return m_bounds;
}

void UIButton::addButtonAction(std::shared_ptr<ButtonAction> &action){
    m_actions.push_back(action);
}

//////*/ for close button only
void UIButton::setWindowToClose(std::string windowID){
    if (m_isCloseButton){
        m_attachedWindow = windowID;
    }

}

bool UIButton::onButtonPress(){
    // if button is a close button, then deactivate everything
    if (m_isCloseButton){
        std::cout << "shownWindowSize: " << m_shownScreens.size() << std::endl;
        //m_shownScreens.erase(m_attachedWindow);
        std::cout << "new shownWindowSize: " << m_shownScreens.size() << std::endl;
        std::cout << "CLOSE WINDOW: " << m_attachedWindow << std::endl;
        return true;


    } else {
        for (auto &action : m_actions){
            action->activate();
        }
        return false;
    }

}


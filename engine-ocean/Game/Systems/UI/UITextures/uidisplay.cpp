#include "uidisplay.h"
#include <set>

UIDisplay::UIDisplay(TextureData tex, glm::vec2 pos, glm::vec2 scale,
                     std::set<std::string>& shownScreens,
                     AspectRatio ratio):
    m_tex(tex),
    m_shownScreens(shownScreens),
    m_aspectRatio(ratio)
{

    m_pos = (pos);
    m_scale = (scale);
    m_tex_aspectRatio = static_cast<float>(m_tex.height)/static_cast<float>(m_tex.width);
    setTransformationMat(pos, scale);
//    std::cout << "tex aspect ratio:" << m_tex_aspectRatio << std::endl;

//    std::cout << "aspect ratio w: " << m_tex.width << std::endl;
//    std::cout << "aspect ratio h: " << m_tex.height << std::endl;


    setWindowPos(640, 480);
//    std::cout << "screen image height: " << m_screenImageHeight  << std::endl;
//    std::cout << "screen image width: " << m_screenImageWidth  << std::endl;

}

UIDisplay::~UIDisplay(){
    glDeleteTextures(1, &m_tex.textureID);
}

void UIDisplay::draw(){}

GLuint UIDisplay::getTexID(){
    return m_tex.textureID;
}

glm::vec2 UIDisplay::getPos(){
 return m_pos;
}

glm::vec2 UIDisplay::getScale(){
 return m_scale;
}

AspectRatio UIDisplay::getAspectRatio(){
    return m_aspectRatio;
}

void UIDisplay::setWindowPos(int width, int height){
    float xpos = .5f*(m_pos.x + 1.f)*static_cast<float>(width);
    float ypos = (1.f-.5f*(m_pos.y + 1.f))*static_cast<float>(height);
    m_windowPos = glm::vec2(xpos, ypos);

    // set everything according to window dimensions
    m_toScreenScale = static_cast<float>(height)*m_scale.y;
    m_screenImageHeight = m_toScreenScale;
    m_screenImageWidth = m_toScreenScale * m_tex_aspectRatio;

    float windowRatio = static_cast<float>(height)/static_cast<float>(width);
    m_textureAspect = windowRatio / m_tex_aspectRatio;
}

void UIDisplay::setPos(glm::vec2 pos){
    m_pos = pos;
    setTransformationMat(m_pos, m_scale);
}
void UIDisplay::setScale(glm::vec2 scale){
    m_scale = scale;
    setTransformationMat(m_pos, m_scale);
}

float UIDisplay::getTextureScaleAspect(){
    return m_textureAspect;
}

glm::vec2 UIDisplay::getWindowPos(){
    return m_windowPos;
}

int UIDisplay::getHeight(){
    return m_tex.height;
}
int UIDisplay::getWidth(){
    return m_tex.width;
}

int UIDisplay::getScreenHeight(){
    return m_screenImageHeight;
}
int UIDisplay::getScreenWidth(){
    return m_screenImageWidth;
}

glm::mat4 UIDisplay::getScaleMatrix(glm::vec2 scale) {
    glm::mat4 M = glm::mat4(1.f);
    M[0][0] = scale.x;//* (m_screenImageHeight/m_screenImageWidth);
    M[1][1] = scale.y; //* (m_tex_aspectRatio);
    M[2][2] = 1.f;
    return M;
}

void UIDisplay::setTransformationMat(glm::vec2 translation, glm::vec2 scale){
    glm::mat4 transMat = glm::mat4(1.f);
    transMat[3] = glm::vec4(translation.x, translation.y, 0.f, 1.f);

    glm::mat4 scaleMat = getScaleMatrix(glm::vec2(scale));

    m_transformationMat = transMat*scaleMat;
}

glm::mat4 UIDisplay::getTransformationMat(){
    return m_transformationMat;
}

// remember that origin is top left corner!!
Bounds2f UIDisplay::getBounds(){
    glm::vec2 halfDimensions = glm::vec2(m_screenImageWidth, m_screenImageHeight)*.5f;
    m_bounds.max = glm::vec2(m_windowPos.x + halfDimensions.x, m_windowPos.y - halfDimensions.y);
    m_bounds.min = glm::vec2(m_windowPos.x - halfDimensions.x, m_windowPos.y + halfDimensions.y);
    return m_bounds;
}

float UIDisplay::getTextureRatio(){
    return m_tex_aspectRatio;
}





#pragma once

// Credit to https://learnopengl.com/In-Practice/Text-Rendering

#include "font.h"
#include <string>
#include <memory>

enum class AnchorPoint{
    TopLeft, TopCenter, TopRight
};

class TextRenderer{
public:
    void initialize();
    void renderUIText(std::shared_ptr<Font> font, std::string text, glm::vec2 anchorPosition, AnchorPoint anchorPoint, float textBoxWidth, float fontSize, float lineSpacing, glm::vec3 textColor);

private:
    GLuint m_vao;
    GLuint m_vbo;
};
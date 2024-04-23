#include "textrenderer.h"
#include "debug.h"

void TextRenderer::initialize(){
    glGenVertexArrays(1, &m_vao);
    glGenBuffers(1, &m_vbo);
    glBindVertexArray(m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);   
}

void TextRenderer::renderUIText(std::shared_ptr<Font> font, std::string text, glm::vec2 anchorPosition, AnchorPoint anchorPoint, float textBoxWidth, float fontSize, float lineSpacing, glm::vec3 textColor){
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(m_vao);

    float x = 0;
    float y = 0;

    float maxY = font->getCharacter(73).size.y * fontSize;
    float minY = - (font->getCharacter(106).size.y - font->getCharacter(106).bearing.y) * fontSize;

    switch(anchorPoint){
        case AnchorPoint::TopLeft:
            x = anchorPosition.x;
            y = anchorPosition.y - maxY;
            break;
        case AnchorPoint::TopCenter:
            x = anchorPosition.x - textBoxWidth * 0.5f;
            y = anchorPosition.y - maxY;
            break;
        case AnchorPoint::TopRight:
            x = anchorPosition.x - textBoxWidth;
            y = anchorPosition.y - maxY;
            break;
    }

    float initX = x;
    float initY = y;

    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++)
    {
        Character ch = font->getCharacter(*c);

        float xpos = x + ch.bearing.x * fontSize;
        float ypos = y - (ch.size.y - ch.bearing.y) * fontSize;

        float w = ch.size.x * fontSize;
        float h = ch.size.y * fontSize;

        if(xpos + w > initX + textBoxWidth){
            x = initX;
            y -= (1.f + lineSpacing) * (maxY - minY);

            xpos = x + ch.bearing.x * fontSize;
            ypos = y - (ch.size.y - ch.bearing.y) * fontSize;
        }


        // update VBO for each character
        float vertices[6][4] = {
            { xpos,     ypos + h,   0.0f, 0.0f },            
            { xpos,     ypos,       0.0f, 1.0f },
            { xpos + w, ypos,       1.0f, 1.0f },

            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos + w, ypos,       1.0f, 1.0f },
            { xpos + w, ypos + h,   1.0f, 0.0f }           
        };
        // render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.textureID);
        // update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); 
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        // render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
        x += (ch.advance >> 6) * fontSize; // bitshift by 6 to get value in pixels (2^6 = 64)
    }
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
}
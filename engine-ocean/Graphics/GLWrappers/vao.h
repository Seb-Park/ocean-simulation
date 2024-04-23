#pragma once

#include "vbo.h"
#include <memory>

enum class VAOAttrib{
    POS = 1,
    NORM = 2,
    UV = 4,
    COLOR = 8
};

inline VAOAttrib operator|(VAOAttrib lhs, VAOAttrib rhs){
    using VAOAttribType = std::underlying_type<VAOAttrib>::type;
    return VAOAttrib(static_cast<VAOAttribType>(lhs) | static_cast<VAOAttribType>(rhs));
}

inline bool operator&(VAOAttrib lhs, VAOAttrib rhs){
    using VAOAttribType = std::underlying_type<VAOAttrib>::type;
    return static_cast<bool>(VAOAttrib(static_cast<VAOAttribType>(lhs) & static_cast<VAOAttribType>(rhs)));
}

class VAO
{
public:
    VAO(std::shared_ptr<VBO> vbo, VAOAttrib attribs);
    ~VAO();

    void bind();
    void unbind();
    
    void draw();

private:
    void addAttribute(GLuint attrib_index, GLint attrib_size);

    std::shared_ptr<VBO> m_vbo; 
    GLuint m_handle;
    GLint m_curr_offset;
    GLuint m_vert_size;

};
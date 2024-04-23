#include "vbo.h"
#include <iostream>

VBO::VBO(std::vector<float> data):
    m_length(data.size())
{
    glGenBuffers(1, &m_handle);
    bind();
    glBufferData(GL_ARRAY_BUFFER, m_length*sizeof(float), data.data(), GL_STATIC_DRAW);
    unbind();
}

VBO::~VBO(){
    glDeleteBuffers(1, &m_handle);
}

void VBO::bind(){
    glBindBuffer(GL_ARRAY_BUFFER, m_handle);
}

void VBO::unbind(){
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

GLsizei VBO::getLength(){
    return m_length;
}
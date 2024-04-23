#include "vao.h"
#include <iostream>

VAO::VAO(std::shared_ptr<VBO> vbo, VAOAttrib attribs):
    m_vbo(vbo),
    m_curr_offset(0),
    m_vert_size(0)
{
    glGenVertexArrays(1, &m_handle);
    bind();

    if(attribs & VAOAttrib::POS){
        m_vert_size += 3;
    }
    if(attribs & VAOAttrib::NORM){
        m_vert_size += 3;
    }
    if(attribs & VAOAttrib::UV){
        m_vert_size += 2;
    }
    if(attribs & VAOAttrib::COLOR){
        m_vert_size += 3;
    }
    

    //Attach layout to vbo depending on which of the above attributes it has
    m_vbo->bind();

    if(attribs & VAOAttrib::POS){
        addAttribute(0, 3);
    }
    if(attribs & VAOAttrib::NORM){
        addAttribute(1, 3);
    }
    if(attribs & VAOAttrib::UV){
        addAttribute(2, 2);
    }
    if(attribs & VAOAttrib::COLOR){
        addAttribute(3, 3);
    }

    m_vbo->unbind();
    unbind();
}

VAO::~VAO(){
    glDeleteVertexArrays(1, &m_handle);
}

void VAO::bind(){
    glBindVertexArray(m_handle);
}

void VAO::unbind(){
    glBindVertexArray(0);
}

void VAO::draw(){
    bind();
    glDrawArrays(GL_TRIANGLES, 0, m_vbo->getLength()/m_vert_size);
    unbind();
}

void VAO::addAttribute(GLuint attrib_index, GLint attrib_size){
    glEnableVertexAttribArray(attrib_index);
    glVertexAttribPointer(attrib_index, attrib_size, GL_FLOAT, GL_FALSE, m_vert_size*sizeof(GLfloat), reinterpret_cast<void*>(m_curr_offset*sizeof(GLfloat)));
    m_curr_offset += attrib_size;
}
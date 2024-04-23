#include "btselector.h"
#include <iostream>

BTSelector::BTSelector()
{

}

void BTSelector::addChildren(BTNode* node){
    m_children.push_back(node);
}

Status BTSelector::update(float seconds){
    // update each children until one doesnt fail
    for (auto node : m_children){
        Status result = node->update(seconds);

        // select this one and return its status --> this node is currently running
        if (result != Status::FAIL){
            m_selected_node = node;
            return result;
        }
    }

    // otherwise if all children fail, then fail this selector
    return Status::FAIL;
}

void BTSelector::reset(){}

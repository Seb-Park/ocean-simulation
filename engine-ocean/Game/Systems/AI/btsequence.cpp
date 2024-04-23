#include "btsequence.h"
#include <vector>

BTSequence::BTSequence()
{

}

void BTSequence::addChildren(BTNode *node){
    m_sequence.push_back(node);
}

Status BTSequence::update(float seconds){

    for (auto node : m_sequence){
        if (node->update(seconds) == Status::FAIL){
            return Status::RUNNING;
        }
        // if come across any node that fails
        if (node->update(seconds) == Status::FAIL){
            return Status::FAIL;
        }
    }

    // if no node is fail or running, sequence is completed
    return Status::SUCCESS;
}

void BTSequence::reset(){}


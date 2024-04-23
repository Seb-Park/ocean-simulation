#ifndef BTSEQUENCE_H
#define BTSEQUENCE_H
#include "btnode.h"
#include <vector>


class BTSequence : public BTNode
{
public:
    BTSequence();
    Status update(float seconds) override;
    void reset() override;
    void addChildren(BTNode* node) override;


private:
    std::vector<BTNode *> m_sequence;
};

#endif // BTSEQUENCE_H

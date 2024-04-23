#ifndef BTSELECTOR_H
#define BTSELECTOR_H
#include "Game/Systems/AI/btnode.h"
#include <vector>


class BTSelector : public BTNode
{
public:
    BTSelector();
    Status update(float seconds) override;
    void reset() override;
    void addChildren(BTNode *node) override;


private:
    std::vector<BTNode *> m_children;
    BTNode *m_selected_node;
};

#endif // BTSELECTOR_H

#ifndef BTCONDITION_H
#define BTCONDITION_H
#include "Game/Systems/AI/btnode.h"


class BTCondition : public BTNode
{
public:
    BTCondition();
    virtual Status update(float seconds) = 0;
    virtual void reset() = 0;


private:
    bool m_condition;
};

#endif // BTCONDITION_H

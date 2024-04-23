#ifndef BTACTION_H
#define BTACTION_H
#include "Game/Systems/AI/btnode.h"

class BTAction : public BTNode
{
public:
    BTAction();
    virtual Status update(float seconds) = 0;
    virtual void reset() = 0;

private:
    Status m_status;

};

#endif // BTACTION_H

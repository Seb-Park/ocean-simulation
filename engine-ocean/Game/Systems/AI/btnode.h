#ifndef BTNODE_H
#define BTNODE_H


enum Status{
    SUCCESS, FAIL, RUNNING
};

class BTNode
{
public:
    BTNode();
    virtual ~BTNode() = default;
    virtual Status update(float seconds) = 0;
    virtual void reset() = 0;
    virtual void addChildren(BTNode *node) = 0;
};
#endif // BTNODE_H

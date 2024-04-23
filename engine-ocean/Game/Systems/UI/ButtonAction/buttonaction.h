#ifndef BUTTONACTION_H
#define BUTTONACTION_H


class ButtonAction
{
public:
    ButtonAction();
    virtual void activate() = 0;
    virtual void deactivate() = 0;
};

#endif // BUTTONACTION_H

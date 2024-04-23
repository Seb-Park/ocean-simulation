#ifndef SHOWWINDOWACTION_H
#define SHOWWINDOWACTION_H
#include "Game/Systems/UI/uisystem.h"
#include <map>
#include "buttonaction.h"

class ShowWindowAction : public ButtonAction
{
public:
    ShowWindowAction(std::map<std::string, std::shared_ptr<UIScreen>>& screens,
                     std::set<std::string>& shownScreens,
                     const std::string screenName);
    void activate() override;
    void deactivate() override;
private:
    std::set<std::string>& m_shownScreens;
    std::string m_screenName;
    std::map<std::string, std::shared_ptr<UIScreen>>& m_screens;
};

#endif // SHOWWINDOWACTION_H

#include "showwindowaction.h"
#include "Game/Systems/UI/uisystem.h"
#include <map>

ShowWindowAction::ShowWindowAction(std::map<std::string, std::shared_ptr<UIScreen>>& all_screens,
                                   std::set<std::string>& shownScreens,
                                   const std::string screenName):
    m_screens(all_screens),
    m_shownScreens(shownScreens)
{
    m_screenName = screenName;
}



void ShowWindowAction::activate(){
    std::cout << "activated window show!!!" << std::endl;

    // add screen to be rendered, and also set it be the only one active
    //m_screens[m_screenName] = m_screen;
//    for (auto &screen : m_screens){
//        screen.second->isActive = false;
//    }

//    m_screens[m_screenName]->isActive = true;
    m_shownScreens.insert(m_screenName);
    m_screens[m_screenName]->isActive = true;

}

void ShowWindowAction::deactivate(){
    m_shownScreens.erase(m_screenName);
}

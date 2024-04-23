#ifndef AIBEHAVIORCOMPONENT_H
#define AIBEHAVIORCOMPONENT_H


#include "Game/Components/Component.h"
#include "Game/Systems/AI/btselector.h"
#include "Game/Systems/aisystem.h"
#include <map>
#include <string>
class AIBehaviorComponent : public Component
{
public:
    AIBehaviorComponent(std::string entity_id,
                        std::map<std::string, BlackboardData>& global_blackboard);
    ~AIBehaviorComponent();
    void update(float seconds);

private:
    void makeBehaviorTree();

    std::string m_entity_id;
    std::map<std::string, BlackboardData>& m_global_blackboard;
    BTNode *m_root = 0;
    Status m_status = Status::SUCCESS;
};

#endif // AIBEHAVIORCOMPONENT_H

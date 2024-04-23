#ifndef PROXIMITYCONDITION_H
#define PROXIMITYCONDITION_H
#include "Game/Components/TransformComponent.h"
#include "btcondition.h"
#include <memory>


class ProximityCondition : public BTNode
{
public:
    ProximityCondition(std::string entity_id,
                       std::map<std::string, BlackboardData>& global_blackboard,
                       float proximity);
    Status update(float seconds) override;
    void reset() override;
    void addChildren(BTNode *node) override;


private:
    bool checkProximity();
    float m_proximity;
    std::string m_entity_id;
    std::map<std::string, BlackboardData>& m_global_blackboard;

};

#endif // PROXIMITYCONDITION_H

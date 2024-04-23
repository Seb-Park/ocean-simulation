#ifndef WALKACTION_H
#define WALKACTION_H
#include "Game/Components/TransformComponent.h"
#include "Game/GameObjects/GameObject.h"
#include "Game/Systems/AI/btnode.h"
#include "glm/fwd.hpp"
#include <memory>
#include "btaction.h"

class WalkAction : public BTNode
{
public:
    WalkAction(std::string entity_id, std::map<std::string, BlackboardData>& global_blackboard);
    Status update(float seconds) override;
    void reset() override;
    void setPath(glm::vec3 entity_pos);
    void addChildren(BTNode *node) override;



private:
    std::vector<glm::vec3> m_path;
    glm::vec3 m_destination;
    std::map<std::string, BlackboardData>& m_global_blackboard;
    std::string m_entity_id;
};


#endif // WALKACTION_H

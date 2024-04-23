#include "aibehaviorcomponent.h"
#include "Game/Systems/AI/Actions/walkaction.h"
#include "Game/Systems/AI/Conditions/proximitycondition.h"
#include "Game/Systems/AI/btselector.h"
#include "Game/Systems/AI/btsequence.h"

AIBehaviorComponent::AIBehaviorComponent(std::string entity_id,
                                         std::map<std::string, BlackboardData>& global_blackboard):
    m_global_blackboard(global_blackboard)
{
    m_entity_id = entity_id;
    makeBehaviorTree();

}

void AIBehaviorComponent::makeBehaviorTree(){
    // leaves
    //std::unique_ptr<BTNode> walk = std::make_unique<WalkAction>(m_entity_id, m_global_blackboard);
    BTNode *proximCond = new ProximityCondition(m_entity_id, m_global_blackboard, 20.f);
    BTNode *walk = new WalkAction(m_entity_id, m_global_blackboard);


//    // pathfind sequence
    BTNode *pathfindSeq = new BTSequence;
    pathfindSeq->addChildren(proximCond);
    pathfindSeq->addChildren(walk);

//    // idle sequence

//    // root
    m_root = new BTSelector;
    m_root->addChildren(pathfindSeq);
}

// how might i be able to generalize the creation of the tree?
void AIBehaviorComponent::update(float seconds){
    // update root, which updates all its children
    //std::cout << "---------in ai system" << std::endl;

    m_status = m_root->update(seconds);
}

AIBehaviorComponent::~AIBehaviorComponent(){
    delete m_root;
}


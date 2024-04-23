#include "ellipsoidtrianglecollisionsystem.h"
#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/CollisionComponents/boundingellipsoid.h"
#include "Game/Components/DrawComponent.h"
#include "Game/Components/TransformComponent.h"
#include "Game/GameObjects/GameObject.h"
#include <memory>

EllipsoidTriangleCollisionSystem::EllipsoidTriangleCollisionSystem(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects) :
    m_rigid_gameobjects(rigid_gameobjects),
    m_bvh_system(std::make_unique<BVHTree>(rigid_gameobjects))
{
}


Triangle EllipsoidTriangleCollisionSystem::convertTriangleToSphereSpace(const glm::vec3 &orig_R, const Triangle &orig_triangle){
    Triangle transformed_triangle;

    // divide vertex pos by Rs
    transformed_triangle.vertexA = orig_triangle.vertexA / orig_R;
    transformed_triangle.vertexB = orig_triangle.vertexB / orig_R;
    transformed_triangle.vertexC = orig_triangle.vertexC / orig_R;

    transformed_triangle.edge1 = transformed_triangle.vertexB - transformed_triangle.vertexA;
    transformed_triangle.edge2 = transformed_triangle.vertexC - transformed_triangle.vertexA;

    // does the normal stay the same across transformations?
    transformed_triangle.normal = glm::normalize(glm::cross(transformed_triangle.edge1, transformed_triangle.edge2));

    return transformed_triangle;
}


glm::vec3 EllipsoidTriangleCollisionSystem::getIntersectionPt(const glm::vec3 &ray_p, const glm::vec3 &ray_d, const float t){
    return ray_p + t*ray_d;
}

float EllipsoidTriangleCollisionSystem::solveQuadratic(float a, float b, float c){
    float discriminant = pow(b,2) - 4*a*c;

    // check for imaginary numbers
    if (discriminant >= 0){
        float pos_quad = (-b + sqrt(discriminant)) / (2*a);
        float neg_quad = (-b - sqrt(discriminant)) / (2*a);

        // if the mininum is a neg number, then return the larger of the two
        if (fmin(pos_quad, neg_quad) < ZERO){
            return fmax(pos_quad, neg_quad);
        } else {
            return fmin(pos_quad, neg_quad);
        }
    }
    return -1;
}


float EllipsoidTriangleCollisionSystem::raycastSphere(const glm::vec3 &ray_p, const glm::vec3 &ray_d, const glm::vec3 posA, const glm::vec3 posB){
    // shift ray and ellipsoid to center first in order to use regular implicit sphere equation
    glm::vec3 p = ray_p - posA;
    glm::vec3 d = ray_d - posA;

    float a = (d.x*d.x) + (d.y*d.y) + (d.z*d.z);
    float b = (2.f*p.x*d.x) + (2.f*p.y*d.y) + (2.f*p.z*d.z);
    float c = (p.x*p.x) + (p.y*p.y) + (p.z*p.z) - 1.f;

    return solveQuadratic(a,b,c);
}

float EllipsoidTriangleCollisionSystem::raycastTriangle(const glm::vec3 &ray_p, const glm::vec3 &ray_d, const Triangle &triangle){
    float t = (glm::dot((triangle.vertexA - ray_p), triangle.normal))/(glm::dot(ray_d, triangle.normal));
    return t;
}

glm::vec3 EllipsoidTriangleCollisionSystem::getTriangleEdge(const glm::vec3 v1, const glm::vec3 v2){
    return v2 - v1;
}

bool EllipsoidTriangleCollisionSystem::isSameDirection(const glm::vec3 edge1, const glm::vec3 edge2, const glm::vec3 n){
    if (glm::dot(glm::cross(edge1, edge2), n) > 0){
        return true;
    }
    return false;
}

void EllipsoidTriangleCollisionSystem::makeCollisionData(const float t, const glm::vec3 triangle_n,
                                                         const glm::vec3 intersection_pos,
                                                         const glm::vec3 point_of_contact,
                                                         const bool isVertex){
    CollisionData data;
    data.t = t;
    data.triangle_n = triangle_n;
    data.intersection_pos = intersection_pos;
    data.point_of_contact = point_of_contact;
    data.isVertex = isVertex;
    m_collisionData.push_back(data);
}


bool EllipsoidTriangleCollisionSystem::detectSphere_with_Interior(const Triangle &triangle, const glm::vec3 posA, const glm::vec3 posB){
    // calculate point on sphere closest to the plane and cast a ray
    glm::vec3 ray_p = posA-triangle.normal; // closest point on unit sphere
    glm::vec3 ray_d = (posB-posA);

    glm::vec3 n = triangle.normal;
    glm::vec3 v = triangle.vertexA;
    float plane_const = -1 * (n.x*v.x + n.y*v.y + n.z*v.z);

    // pdf equation
    float t = (1-(glm::dot(n, posA) + plane_const)) / (glm::dot(n, (posB-posA)));
    if (t < ZERO || t > ONE){
        return false;
    }

    // the point P on the sphere closest to plane surface
    glm::vec3 P_edge = getIntersectionPt(ray_p, ray_d, t);


    // determine if point is on triangle plane
        // determine if point is also within triangle bounds
        glm::vec3 AB = getTriangleEdge(triangle.vertexA, triangle.vertexB);
        glm::vec3 BC = getTriangleEdge(triangle.vertexB, triangle.vertexC);
        glm::vec3 CA = getTriangleEdge(triangle.vertexC, triangle.vertexA);

        glm::vec3 AP = getTriangleEdge(triangle.vertexA, P_edge);
        glm::vec3 BP = getTriangleEdge(triangle.vertexB, P_edge);
        glm::vec3 CP = getTriangleEdge(triangle.vertexC, P_edge);

        if (isSameDirection(AB, AP, triangle.normal) &&
                isSameDirection(BC, BP, triangle.normal) &&
                isSameDirection(CA, CP, triangle.normal)){

            // the actual center of sphere
            glm::vec3 P_centered = getIntersectionPt(posA, posB-posA, t);
            makeCollisionData(t, triangle.normal, P_centered, P_edge, false);
            return true;
        }
    return false;
}

// do this for each vertex pair AB, BC, CA
bool EllipsoidTriangleCollisionSystem::detectSphere_with_Edge(const glm::vec3 &v1, const glm::vec3 &v2,
                                                              const glm::vec3 &triangle_normal,
                                                              const glm::vec3 &posA, const glm::vec3 &posB){
    glm::vec3 C = v1;
    glm::vec D = v2;

    float a = pow(glm::length(glm::cross(posB-posA, D-C)),2);
    float b = glm::dot(2.f*glm::cross(posB-posA, D-C), glm::cross(posA-C, D-C));
    float c = pow(glm::length(glm::cross(posA-C, D-C)),2) - pow(glm::length(D-C),2);

    float t = solveQuadratic(a,b,c);

    // discard a negative t --> indicates no collision with edge
    if (t < ZERO || t > ONE){
        return false;
    }

    glm::vec3 P = getIntersectionPt(posA, posB-posA, t);

    glm::vec3 AP = P-C;
    glm::vec3 AB = D-C;
    glm::vec3 normal = AP - (glm::dot(AB,AP)/glm::dot(AB,AB))*AB;

    // determine if intersection pt is between C and D
    if ( (0 < glm::dot(P-C,D-C)) &&  (glm::dot(P-C,D-C) < pow(glm::length(D-C),2)) ){
        makeCollisionData(t, normal, P, P, false);
        return true;
    }

    return false;
}

// needs to be done with each vertex
bool EllipsoidTriangleCollisionSystem::detectSphere_with_Vertex(const glm::vec3 &v, const glm::vec3 triangle_normal, const glm::vec3 &posA, const glm::vec3 &posB){
    glm::vec3 ray_p = v;
    //glm::vec3 ray_d = v-(posB-posA);

    // convert to origin space:
    glm::vec3 vel = posB-posA;
    glm::vec3 v_origined = v - posA;

    float a = glm::dot(vel,vel);
    float b = 2*(glm::dot(vel, -v_origined));
    float c = pow(glm::length(v_origined),2) - 1.f;

    float t = solveQuadratic(a,b,c);
    if (t < ZERO || t > ONE){
        return false;
    }

    // where the center is when the sphere hits vertex
    glm::vec3 center_pt_origined = t*vel;
    // convert back to non-origin space
    glm::vec3 center_pt = center_pt_origined + posA;
    glm::vec3 normal = center_pt - v;
    makeCollisionData(t, normal, center_pt, v, true);
     return true;
}

CollisionData EllipsoidTriangleCollisionSystem::getNearestCollision(const std::vector<CollisionData> &collisionData, const Triangle &sphere_tri){
    CollisionData nearestCollision;
    nearestCollision.t = INFINITY;

    for (const CollisionData &data : collisionData){
        if (data.t < nearestCollision.t){
            // set nearest collision to datapoint
            nearestCollision = data;
        }
    }
    return nearestCollision;
}

// SPHERE SPACE
CollisionData EllipsoidTriangleCollisionSystem::detectTriangleCollision(const Triangle &sphere_tri, const glm::vec3 &posA, const glm::vec3 &posB){

    m_collisionData.clear();
    // if interior is true, no need to check rest. get the only collision data
    if (detectSphere_with_Interior(sphere_tri, posA, posB)){
        return m_collisionData[0];
    }

    // otherwise check rest:

    // detect edge with all 3 edges
    // | operator executes regardless of if left hand side is true or not
    if (
    detectSphere_with_Edge(sphere_tri.vertexA, sphere_tri.vertexB, sphere_tri.normal, posA, posB) |
    detectSphere_with_Edge(sphere_tri.vertexB, sphere_tri.vertexC, sphere_tri.normal, posA, posB) |
    detectSphere_with_Edge(sphere_tri.vertexC, sphere_tri.vertexA, sphere_tri.normal, posA, posB) ){
        return getNearestCollision(m_collisionData, sphere_tri);
    }

    // detect vertex with all 3 vertices
    detectSphere_with_Vertex(sphere_tri.vertexA, sphere_tri.normal, posA, posB);
    detectSphere_with_Vertex(sphere_tri.vertexB, sphere_tri.normal, posA, posB);
    detectSphere_with_Vertex(sphere_tri.vertexC, sphere_tri.normal, posA, posB);

    return getNearestCollision(m_collisionData, sphere_tri);
}

CollisionData EllipsoidTriangleCollisionSystem::revertDataToOriginalSpace(const CollisionData &sphere_datum, const glm::vec3 &ellip_R, const glm::vec3 orig_triangle_n){
    CollisionData reverted_datum = sphere_datum;
    reverted_datum.intersection_pos = reverted_datum.intersection_pos*ellip_R;
    reverted_datum.point_of_contact = reverted_datum.point_of_contact*ellip_R;

    return reverted_datum;
}

void EllipsoidTriangleCollisionSystem::setEllipsoidDimensions(const glm::vec3 &dimensions){
    m_ellipsoid_dimensions = dimensions;
}

CollisionData EllipsoidTriangleCollisionSystem::collideWithWorld(const glm::vec3 &posA, const glm::vec3 &posB){
    // makes ellipsoids
    Ellipsoid ellipsoidA;
    ellipsoidA.R = m_ellipsoid_dimensions;
    ellipsoidA.center_pos = posA;


    Ellipsoid ellipsoidB;
    ellipsoidB.R = m_ellipsoid_dimensions;
    ellipsoidB.center_pos = posB;

    // convert ellipsoid centers to sphere space
    glm::vec3 sphere_posA = ellipsoidA.center_pos/m_ellipsoid_dimensions;
    glm::vec3 sphere_posB = ellipsoidB.center_pos/m_ellipsoid_dimensions;

    CollisionData nearestCollision;
    nearestCollision.t = INFINITY;


    std::vector<Triangle> candidates = m_bvh_system->getBVHDetectedCollisions(posA, posB, ellipsoidA.R);
   // std::cout << "candidate size: " << candidates.size() << std::endl;


//    for (auto &rigid_go : candidates){
//        //std::cout << "candidate shape name: " << rigid_go->getComponent<DrawComponent>()->getShapeName() << std::endl;
//        if (rigid_go->hasComponent<CollisionComponent>()){
//            if (rigid_go->getComponent<CollisionComponent>()->hasCollisionShape<BoundingTriangle>()){
//                // for every triangle in game obj, transform to sphere space and collide with ellipsoid (aka sphere)
                for (int i=0; i< candidates.size(); i++){
                    Triangle triangle = candidates[i];
//                    std::cout<< "triangle i: " << i << ", v1: (" << triangle.vertexA.x << ", " << triangle.vertexA.y << ", " << triangle.vertexA.z << "), v2: ("
//                                                                 << triangle.vertexB.x << ", " << triangle.vertexB.y << ", " << triangle.vertexB.z << "), v3: ("
//                                                                 << triangle.vertexC.x << ", " << triangle.vertexC.y << ", " << triangle.vertexC.z << ")" << std::endl;

                    // transform triangle to sphere space and detect collision data (in sphere space) for this specific triangle
                    Triangle sphere_triangle = convertTriangleToSphereSpace(m_ellipsoid_dimensions, triangle);
                    CollisionData sphere_collisionDatum = detectTriangleCollision(sphere_triangle, sphere_posA, sphere_posB);
                    CollisionData ellip_collisionDatum = revertDataToOriginalSpace(sphere_collisionDatum, m_ellipsoid_dimensions, triangle.normal);

                    // keep track of nearest Collision
                    if (ellip_collisionDatum.t < nearestCollision.t && ellip_collisionDatum.t <= ONE && ellip_collisionDatum.t >= ZERO){
                        nearestCollision = ellip_collisionDatum;
                    }
                    // repeat for the other triangles
                }
            //}
        //}
        // repeat for all gameobjects
    //}

    return nearestCollision;
}

glm::vec3 EllipsoidTriangleCollisionSystem::doNudge(glm::vec3 curr_pos, CollisionData collision){
    glm::vec3 nudge = m_scale_normals_mt*collision.triangle_n;

    glm::vec3 pos_nudged =  collision.intersection_pos + nudge * .01f;
    int MAX_NUDGES = 3;
    for (int i=0; i<MAX_NUDGES; i++){
        CollisionData nudge_collision = collideWithWorld(curr_pos, pos_nudged);
        if (nudge_collision.t == INFINITY){
            curr_pos = pos_nudged;
            break;
        } else {
            if (glm::length((m_scale_normals_mt*nudge_collision.triangle_n) - nudge) < .0001 || glm::length((m_scale_normals_mt*nudge_collision.triangle_n) + nudge) < .0001){
                nudge = -m_scale_normals_mt*nudge_collision.triangle_n;
            } else {
                nudge = m_scale_normals_mt*nudge_collision.triangle_n;
            }
            pos_nudged = nudge_collision.intersection_pos + nudge * .01f;
        }
    }
    return curr_pos;

}

glm::mat3 EllipsoidTriangleCollisionSystem::makeConversionMT(glm::vec3 initial_center_pos, glm::vec3 dimensions){
    m_conversion_mt = std::make_shared<ModelTransform>();
    m_conversion_mt->translate(-initial_center_pos);
    m_conversion_mt->scale(dimensions);
    m_scale_normals_mt = glm::transpose((glm::inverse(glm::mat3(m_conversion_mt->getModelMatrix()))));
    m_scale_normals_mt = glm::mat4(1.f);
    return m_scale_normals_mt;
}

// in sphere space
std::pair<std::vector<CollisionData>, glm::vec3> EllipsoidTriangleCollisionSystem::mtvSlide(glm::vec3 initial_pos, glm::vec3 final_pos, glm::vec3 dimensions){
    setEllipsoidDimensions(dimensions);
    glm::mat3 normals_scale = makeConversionMT(initial_pos, dimensions);

    std::vector<CollisionData> collisions;
    glm::vec3 curr_pos = initial_pos;
    glm::vec3 next_pos = final_pos;
    int MAX_TRANSLATIONS = 3;
    for (int i=0; i<MAX_TRANSLATIONS; i++){
        CollisionData c = collideWithWorld(initial_pos, next_pos);

        // if no collisions [0,1] found, then continue on regular path
        if (c.t == INFINITY){
            //std::cout << "--mtv-- " << next_pos.x << "," << next_pos.y << "," << next_pos.z << std::endl;
            return std::pair<std::vector<CollisionData>, glm::vec3>(collisions, next_pos);
        } else {
            // nudge along the collision data
            //std::cout << "COLLIDING......" << std::endl;
            curr_pos = doNudge(curr_pos, c);
            glm::vec3 d = next_pos - curr_pos;
            glm::vec3 d_corrected = (d - glm::dot(d, normals_scale*c.triangle_n)) * (normals_scale*c.triangle_n);
            next_pos = curr_pos + d_corrected;
            collisions.push_back(c);
        }

    }

    return std::pair<std::vector<CollisionData>, glm::vec3>(collisions, curr_pos);
}



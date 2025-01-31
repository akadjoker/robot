#include "RigidBody.hpp"



std::vector<RigidBody *> bodies;

extern void RemoveAllJoints();

void ClearWorld()
{
    for (auto body : bodies)
    {
        delete body;
    }
    RemoveAllJoints();
    bodies.clear();
}

void AddBody(RigidBody *body)
{
    bodies.push_back(body);
}

RigidBody::RigidBody(bool isStatic)
{
    this->id = 0;
    this->isStatic = isStatic;
    this->body = nullptr;

    if (!isStatic)
    {
        body = dBodyCreate(world);
    }
    AddBody(this);
    printf("RigidBody created\n");
}

RigidBody::~RigidBody()
{

    dGeomDestroy(geom);
    if (body)
    {
        dBodyDestroy(body);
    }
    body = nullptr;
    geom = nullptr;
    printf("RigidBody destroyed\n");
}

void RigidBody::setGeomPosition(dReal x, dReal y, dReal z)
{
    if (body)
    {
        dGeomSetPosition(geom, x, y, z);
    }
}

void RigidBody::setBodyPosition(dReal x, dReal y, dReal z)
{
    if (body)
    {
        dBodySetPosition(body, x, y, z);
    }
}

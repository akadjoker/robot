#pragma once
#include "ode/ode.h"
#include "RigidBody.hpp"

class Joint
{
protected:
    dJointID joint;
    dJointType type;
    RigidBody *body1;
    RigidBody *body2;
    long id;

public:
    Joint();
    virtual ~Joint();
    void setType(dJointType type);
    dJointID getJoint();
    void attach(RigidBody *body1, RigidBody *body2);

    void setId(long id) { this->id = id; }
    long getId() { return id; }
};

//  (Dobradiça)
class HingeJoint : public Joint
{
public:
    HingeJoint();
    void setAnchor(dReal x, dReal y, dReal z);
    void setAxis(dReal x, dReal y, dReal z);
    void setLimits(dReal lo, dReal hi);
};

//  Hinge2 Joint (Para Rodas)
class Hinge2Joint : public Joint
{
public:
    Hinge2Joint();
    void setAnchor(dReal x, dReal y, dReal z);

    void setAxes(dReal ax, dReal ay, dReal az, dReal bx, dReal by, dReal bz);

    void setSuspension(dReal ERP, dReal CFM);

    void setAxis1(dReal x, dReal y, dReal z);
    void setAxis2(dReal x, dReal y, dReal z);
};

class FixedJoint : public Joint
{
public:
    FixedJoint();
    void setFixed();
};

class SliderJoint : public Joint
{
public:
    SliderJoint();
    void setAxis(dReal x, dReal y, dReal z);
    void setLimits(dReal lo, dReal hi);
};

// 🔹 Ball Joint (Articulação Esférica)
class BallJoint : public Joint
{
public:
    BallJoint();
    void setAnchor(dReal x, dReal y, dReal z);
};

// 🔹 AMotor Joint (Motor Angular)
class AMotorJoint : public Joint
{
public:
    AMotorJoint();
    void setNumAxes(int num);
    void setAxis(int axisIndex, int relativeTo, dReal x, dReal y, dReal z);
    void setVelocity(int axisIndex, dReal velocity);
    void setFMax(int axisIndex, dReal maxForce);
};

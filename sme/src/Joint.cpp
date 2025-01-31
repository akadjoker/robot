#include "Joint.hpp"
#include <vector>

 std::vector<Joint *> joints;
extern dWorldID world;
extern dSpaceID space;


// dJointID dJointCreateBall (dWorldID, dJointGroupID);
// dJointID dJointCreateHinge (dWorldID, dJointGroupID);
// dJointID dJointCreateSlider (dWorldID, dJointGroupID);
// dJointID dJointCreateContact (dWorldID, dJointGroupID, const dContact *);
// dJointID dJointCreateUniversal (dWorldID, dJointGroupID);
// dJointID dJointCreateHinge2 (dWorldID, dJointGroupID);
// dJointID dJointCreatePR (dWorldID, dJointGroupID);
// dJointID dJointCreatePU (dWorldID, dJointGroupID);
// dJointID dJointCreatePiston (dWorldID, dJointGroupID);
// dJointID dJointCreateFixed (dWorldID, dJointGroupID);
// dJointID dJointCreateAMotor (dWorldID, dJointGroupID);
// dJointID dJointCreateLMotor (dWorldID, dJointGroupID);
// dJointID dJointCreatePlane2D (dWorldID, dJointGroupID);
// dJointID dJointCreateDBall (dWorldID, dJointGroupID);
// dJointID dJointCreateDHinge (dWorldID, dJointGroupID);
// dJointID dJointCreateTransmission (dWorldID, dJointGroupID);

void RemoveAllJoints()
{
    for (auto joint : joints)
     {
         delete joint;
     }
    joints.clear();
}

Joint::Joint() : joint(nullptr), type(dJointTypeNone), body1(nullptr), body2(nullptr), id(0) 
{
    joints.push_back(this);
    printf("Create joint\n");
}

Joint::~Joint()
{
    if (joint)
    {
        dJointDestroy(joint);
    }
 printf("Destroy joint\n");
}

void Joint::setType(dJointType type)
{
    this->type = type;
  
}

dJointID Joint::getJoint()
{
    return joint;
}

void Joint::attach(RigidBody *rb1, RigidBody *rb2)
{
    body1 = rb1;
    body2 = rb2;
    dJointAttach(joint, body1->getBody(), body2->getBody());
}


HingeJoint::HingeJoint()
{
    setType(dJointTypeHinge);
}

void HingeJoint::setAnchor(dReal x, dReal y, dReal z)
{
    dJointSetHingeAnchor(joint, x, y, z);
}

void HingeJoint::setAxis(dReal x, dReal y, dReal z)
{
    dJointSetHingeAxis(joint, x, y, z);
}

void HingeJoint::setLimits(dReal lo, dReal hi)
{
    dJointSetHingeParam(joint, dParamLoStop, lo);
    dJointSetHingeParam(joint, dParamHiStop, hi);
}


Hinge2Joint::Hinge2Joint()
{
    setType(dJointTypeHinge2);
}

void Hinge2Joint::setAnchor(dReal x, dReal y, dReal z)
{
    dJointSetHinge2Anchor(joint, x, y, z);
}


void Hinge2Joint::setAxes(dReal ax, dReal ay, dReal az, dReal bx, dReal by, dReal bz)
{
    dReal a[3] = {ax, ay, az};
    dReal b[3] = {bx, by, bz};

    dJointSetHinge2Axes(joint, a, b);
}

void Hinge2Joint::setSuspension(dReal ERP, dReal CFM)
{
    dJointSetHinge2Param(joint, dParamSuspensionERP, ERP);
    dJointSetHinge2Param(joint, dParamSuspensionCFM, CFM);
}

// void Hinge2Joint::setAxis1(dReal x, dReal y, dReal z)
// {
//     dJointSetHinge2Axis1(joint, x, y, z);
// }

// void Hinge2Joint::setAxis2(dReal x, dReal y, dReal z)
// {
//     dJointSetHinge2Axis2(joint, x, y, z);
// }

FixedJoint::FixedJoint()
{
    setType(dJointTypeFixed);
}

void FixedJoint::setFixed()
{
    dJointSetFixed(joint);
}

SliderJoint::SliderJoint()
{
    setType(dJointTypeSlider);
}

void SliderJoint::setAxis(dReal x, dReal y, dReal z)
{
    dJointSetSliderAxis(joint, x, y, z);
}

void SliderJoint::setLimits(dReal lo, dReal hi)
{
    dJointSetSliderParam(joint, dParamLoStop, lo);
    dJointSetSliderParam(joint, dParamHiStop, hi);
}


BallJoint::BallJoint()
{
    setType(dJointTypeBall);
}

void BallJoint::setAnchor(dReal x, dReal y, dReal z)
{
    dJointSetBallAnchor(joint, x, y, z);
}

// 🔹 Implementação de `AMotorJoint`
AMotorJoint::AMotorJoint()
{
    setType(dJointTypeAMotor);
}

void AMotorJoint::setNumAxes(int num)
{
    dJointSetAMotorNumAxes(joint, num);
}

void AMotorJoint::setAxis(int axisIndex, int relativeTo, dReal x, dReal y, dReal z)
{
    dJointSetAMotorAxis(joint, axisIndex, relativeTo, x, y, z);
}

void AMotorJoint::setVelocity(int axisIndex, dReal velocity)
{
    dJointSetAMotorParam(joint, dParamVel + axisIndex, velocity);
}

void AMotorJoint::setFMax(int axisIndex, dReal maxForce)
{
    dJointSetAMotorParam(joint, dParamFMax + axisIndex, maxForce);
}

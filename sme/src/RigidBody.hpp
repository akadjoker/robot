#pragma once
#include "ode/ode.h"
#include <vector>
extern dWorldID world;
extern dSpaceID space;

void ClearWorld();

class RigidBody
{
protected:
    dBodyID body = nullptr;
    dGeomID geom = nullptr;
    long id;
    bool isStatic;

public:
    RigidBody(bool isStatic);

    virtual ~RigidBody();

    void setGeomPosition(dReal x, dReal y, dReal z);
    void setBodyPosition(dReal x, dReal y, dReal z);

    void setId(long id) { this->id = id; }
    long getId() { return id; }

    void setPosition(dReal x, dReal y, dReal z)
    {
        if (geom)
            dGeomSetPosition(geom, x, y, z); //  estáticos só têm `dGeomSetPosition`
    }

    void setMass(dMass &mass)
    {
        if (body)
        {
            dBodySetMass(body, &mass);
        }
    }

    dBodyID getBody() { return body; }
    dGeomID getGeom() { return geom; }
};

class Cube : public RigidBody
{
public:
    Cube(dReal w, dReal h, dReal d, dReal massValue, bool isStatic)
        : RigidBody(isStatic)
    {
        geom = dCreateBox(space, w, h, d);

        if (!isStatic)
        {
            dMass mass;
            dMassSetBox(&mass, 1.0, w, h, d);
            dMassAdjust(&mass, massValue);
            setMass(mass);
            dGeomSetBody(geom, body);
        }
    }
};

class Sphere : public RigidBody
{
public:
    Sphere(dReal radius, dReal massValue, bool isStatic)
        : RigidBody(isStatic)
    {
        geom = dCreateSphere(space, radius);

        if (!isStatic)
        {
            dMass mass;
            dMassSetSphere(&mass, 1.0, radius);
            dMassAdjust(&mass, massValue);
            setMass(mass);
            dGeomSetBody(geom, body);
        }
    }
};

class Capsule : public RigidBody
{
public:
    Capsule(int direction, dReal radius, dReal length, dReal massValue, bool isStatic)
        : RigidBody(isStatic)
    {
        geom = dCreateCapsule(space, radius, length);

        if (!isStatic)
        {
            dMass mass;
            dMassSetCapsule(&mass, 1, direction, radius, length); // Direção = 3 (eixo Z)
            dMassAdjust(&mass, massValue);
            setMass(mass);
            dGeomSetBody(geom, body);
        }
    }
};

class Cylinder : public RigidBody
{
public:
    Cylinder(int direction, dReal radius, dReal length, dReal massValue, bool isStatic = false)
        : RigidBody(isStatic)
    {
        geom = dCreateCylinder(space, radius, length);

        if (!isStatic)
        {
            dMass mass;
            dMassSetCylinder(&mass, 1, direction, radius, length); // Direção = 3 (eixo Z)
            dMassAdjust(&mass, massValue);
            setMass(mass);
            dGeomSetBody(geom, body);
        }
    }
};

class Vehicle : public RigidBody
{
public:
    std::vector<dBodyID> wheels;
    std::vector<dGeomID> wheelGeoms;
    std::vector<dJointID> wheelJoints;

    Vehicle(dReal mass, dReal width, dReal height, dReal length)
        : RigidBody(false)
    {

        printf("Vehicle created\n");

        dMass m;
        dMassSetBox(&m, 1, width, height, length); // Densidade
        dMassAdjust(&m, mass);                     // Definir massa

        dBodySetMass(body, &m);
        dBodySetAutoDisableFlag(body, 0); // Nunca desativar automaticamente

        geom = dCreateBox(space, width, height, length);
        dGeomSetBody(geom, body);
    }

    void addWheel(float radius, float width, float mass, float x, float y, float z, float ax, float ay, float az)
    {
        dMass m;
        dMassSetCylinder(&m, 1, 3, radius, width); // Direção do eixo = 3 (Z)
        dMassAdjust(&m, mass);

        dQuaternion q;
        dQFromAxisAndAngle(q, 0, 0, 1, M_PI * 0.5);



        dBodyID wheelBody = dBodyCreate(world);
        dBodySetMass(wheelBody, &m);
        dBodySetQuaternion(wheelBody, q);

        dGeomID wheelGeom = dCreateCylinder(space, radius, width);
        dGeomSetBody(wheelGeom, wheelBody);
        dBodySetFiniteRotationMode( wheelBody, 1);
        dBodySetAutoDisableFlag( wheelBody, 0);

        const dReal *cp = dBodyGetPosition(body);
        dReal jx = cp[0] + x;
        dReal jy = cp[1] + y;
        dReal jz = cp[2] + z;

        dBodySetPosition(wheelBody, jx, jy, jz);
        
       

        // Criar hinge2 joint (suspensão + direção)
        dJointID joint = dJointCreateHinge2(world, 0);
        dJointAttach(joint, body, wheelBody);
        dJointSetHinge2Anchor( joint,  jx, jy, jz);

        dReal axis1[] = {0, 1, 0};                                  // Direção de rotação principal (eixo Y)
        dReal axis2[] = { ax, ay, az}; // Rotação da roda

        dJointSetHinge2Axes(joint, axis1, axis2);



        // Guardar referências
        wheels.push_back(wheelBody);
        wheelGeoms.push_back(wheelGeom);
        wheelJoints.push_back(joint);
    }

    void SetParam(int index,int type, float value)
    {
        if (index < wheels.size())
        {
            dJointSetHinge2Param(wheelJoints[index],  type, value);
        }
             
    }

    void setSuspension(int index, dReal ERP, dReal CFM)
    {
        if (index < wheels.size())
        {
            dJointSetHinge2Param(wheelJoints[index], dParamSuspensionERP, ERP);
            dJointSetHinge2Param(wheelJoints[index], dParamSuspensionCFM, CFM);
        }
    }
    void setAcceleration(int index, float accel, float maxAccelForce)
    {
        if (index < wheels.size())
        {
            float target = (fabs(accel) > 0.1) ? maxAccelForce : 0;

            dJointSetHinge2Param(wheelJoints[index], dParamVel2, (dReal)accel);
            dJointSetHinge2Param(wheelJoints[index], dParamFMax2, (dReal)target);
        }
    }

    void setSteering(int index, float steer, float steerFactor)
    {
        if (index < wheels.size())
        {
            dReal v = steer - dJointGetHinge2Angle1(wheelJoints[index]);
            v *= steerFactor;
            dJointSetHinge2Param(wheelJoints[index], dParamVel, (dReal)v);
        }
    }

    const dReal *getLinearVelocity() const
    {
        return dBodyGetLinearVel(body);
    }

    // Obter velocidade angular do veículo
    const dReal *getAngularVelocity() const
    {
        return dBodyGetAngularVel(body);
    }

    float getWheelSpeed(int index) const
    {
        if (index < wheels.size())
        {
            return (float)dJointGetHinge2Angle2Rate(wheelJoints[index]);
        }
        return 0;
    }

    float getSteeringAngle(int index) const
    {
        if (index < wheels.size())
        {
            return (float)dJointGetHinge2Angle1(wheelJoints[index]);
        }
        return 0;
    }

    void applyBrake(int index, dReal brakingForce)
    {
        if (index < wheels.size())
        {
            dJointSetHinge2Param(wheelJoints[index], dParamVel2, 0);
            dJointSetHinge2Param(wheelJoints[index], dParamFMax2, brakingForce);
        }
    }

    void applyBrakesAll(dReal brakingForce)
    {
        for (size_t i = 0; i < wheels.size(); i++)
        {
            applyBrake(i, brakingForce);
        }
    }

    // Verificar se o veículo está de rodas para cima
    bool isUpsideDown() const
    {
        const dReal *rot = dBodyGetRotation(body);
        // Vetor Y do chassis (rot[1], rot[5], rot[9])
        // Se o Y local estiver apontando para baixo, o carro está capotado
        return rot[5] < 0;
    }

    void reset(float x, float y, float z,float angle, float ax, float ay, float az)
    {

        dBodySetPosition(body, x, y, z);

        dMatrix3 R;
        dRSetIdentity(R);
        dRFromAxisAndAngle(R, ax, ay, az, angle);
        dBodySetRotation(body, R);

        dBodySetLinearVel(body, 0, 0, 0);
        dBodySetAngularVel(body, 0, 0, 0);


        for (size_t i = 0; i < wheels.size(); i++)
        {
            dReal jointPos[3];
            dJointGetHinge2Anchor(wheelJoints[i], jointPos);
            dBodySetPosition(wheels[i], jointPos[0], jointPos[1], jointPos[2]);
            dBodySetLinearVel(wheels[i], 0, 0, 0);
            dBodySetAngularVel(wheels[i], 0, 0, 0);
        }
    }

    ~Vehicle()
    {
        printf("Destroy vehicle\n");
        for (size_t i = 0; i < wheels.size(); ++i)
        {
            dBodyDestroy(wheels[i]);
            dGeomDestroy(wheelGeoms[i]);
            dJointDestroy(wheelJoints[i]);
        }
    }
};

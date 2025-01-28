#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
#include "ode/ode.h"

void getTransformationMatrixFromQuaternion(const dReal *p, const dQuaternion q, Matrix *mat)
{

    Quaternion rotation;
    rotation.w = (float)q[0];
    rotation.x = (float)q[1];
    rotation.y = (float)q[2];
    rotation.z = (float)q[3];

    *mat = QuaternionToMatrix(rotation);

    mat->m12 = p[0];
    mat->m13 = p[1];
    mat->m14 = p[2];
    mat->m15 = 1;
}

void RenderLine(const dVector3 p1, const dVector3 p2, Color c)
{
    float p1f[3], p2f[3];
    for (int i = 0; i < 3; i++)
    {
        p1f[i] = (float)p1[i];
        p2f[i] = (float)p2[i];
    }
    DrawLine3D((Vector3){p1f[0], p1f[1], p1f[2]}, (Vector3){p2f[0], p2f[1], p2f[2]}, c);
}
void DrawDetailedCylinder(float length, float radius, Color color)
{
    const int segments = 24; //  segmentos (divisível por 4)
    float halfLength = length * 0.5f;
    float angleStep = (2.0f * PI) / segments;

    rlBegin(RL_TRIANGLES);
    for (int i = 0; i < segments; i++)
    {
        float angle = angleStep * i;
        float nextAngle = angleStep * (i + 1);

        float x1 = cos(angle) * radius;
        float y1 = sin(angle) * radius;
        float x2 = cos(nextAngle) * radius;
        float y2 = sin(nextAngle) * radius;

        rlColor4ub(color.r, color.g, color.b, color.a);
        rlVertex3f(x1, y1, -halfLength);
        rlVertex3f(x2, y2, -halfLength);
        rlVertex3f(x1, y1, halfLength);

        rlVertex3f(x2, y2, -halfLength);
        rlVertex3f(x2, y2, halfLength);
        rlVertex3f(x1, y1, halfLength);

        rlColor4ub(color.r * 0.9, color.g * 0.9, color.b * 0.9, color.a);
        rlVertex3f(0, 0, halfLength);
        rlVertex3f(x1, y1, halfLength);
        rlVertex3f(x2, y2, halfLength);

        rlColor4ub(color.r * 0.8, color.g * 0.8, color.b * 0.8, color.a);
        rlVertex3f(0, 0, -halfLength);
        rlVertex3f(x2, y2, -halfLength);
        rlVertex3f(x1, y1, -halfLength);
    }
    rlEnd();

    rlBegin(RL_LINES);
    for (int i = 0; i < segments; i++)
    {
        float angle = angleStep * i;
        float nextAngle = angleStep * (i + 1);

        float x1 = cos(angle) * radius;
        float y1 = sin(angle) * radius;
        float x2 = cos(nextAngle) * radius;
        float y2 = sin(nextAngle) * radius;

        rlColor4ub(255, 255, 255, 255);
        rlVertex3f(x1, y1, -halfLength);
        rlVertex3f(x2, y2, -halfLength);
        rlVertex3f(x1, y1, halfLength);
    }
    rlEnd();
}

void DrawDetailedCapsule(float length, float radius, Color color)
{
    const int segments = 24;
    const int rings = 12; //  semiesferas
    float halfLength = length * 0.5f;

    DrawDetailedCylinder(length, radius, color);

    rlBegin(RL_TRIANGLES);

    for (int i = 0; i < segments; i++)
    {
        for (int j = 0; j < rings; j++)
        {
            float phi1 = PI * 0.5f * j / rings;
            float phi2 = PI * 0.5f * (j + 1) / rings;
            float theta1 = 2.0f * PI * i / segments;
            float theta2 = 2.0f * PI * (i + 1) / segments;

            float x1 = radius * cos(phi1) * cos(theta1);
            float y1 = radius * cos(phi1) * sin(theta1);
            float z1 = radius * sin(phi1) + halfLength;

            float x2 = radius * cos(phi1) * cos(theta2);
            float y2 = radius * cos(phi1) * sin(theta2);
            float z2 = radius * sin(phi1) + halfLength;

            float x3 = radius * cos(phi2) * cos(theta2);
            float y3 = radius * cos(phi2) * sin(theta2);
            float z3 = radius * sin(phi2) + halfLength;

            float x4 = radius * cos(phi2) * cos(theta1);
            float y4 = radius * cos(phi2) * sin(theta1);
            float z4 = radius * sin(phi2) + halfLength;

            rlColor4ub(color.r * 0.9, color.g * 0.9, color.b * 0.9, color.a);
            rlVertex3f(x1, y1, z1);
            rlVertex3f(x2, y2, z2);
            rlVertex3f(x3, y3, z3);

            rlVertex3f(x1, y1, z1);
            rlVertex3f(x3, y3, z3);
            rlVertex3f(x4, y4, z4);
        }
    }

    for (int i = 0; i < segments; i++)
    {
        for (int j = 0; j < rings; j++)
        {
            float phi1 = PI * 0.5f * j / rings;
            float phi2 = PI * 0.5f * (j + 1) / rings;
            float theta1 = 2.0f * PI * i / segments;
            float theta2 = 2.0f * PI * (i + 1) / segments;

            float x1 = radius * cos(phi1) * cos(theta1);
            float y1 = radius * cos(phi1) * sin(theta1);
            float z1 = -(radius * sin(phi1) + halfLength);

            float x2 = radius * cos(phi1) * cos(theta2);
            float y2 = radius * cos(phi1) * sin(theta2);
            float z2 = -(radius * sin(phi1) + halfLength);

            float x3 = radius * cos(phi2) * cos(theta2);
            float y3 = radius * cos(phi2) * sin(theta2);
            float z3 = -(radius * sin(phi2) + halfLength);

            float x4 = radius * cos(phi2) * cos(theta1);
            float y4 = radius * cos(phi2) * sin(theta1);
            float z4 = -(radius * sin(phi2) + halfLength);

            rlColor4ub(color.r * 0.8, color.g * 0.8, color.b * 0.8, color.a);
            rlVertex3f(x1, y1, z1);
            rlVertex3f(x3, y3, z3);
            rlVertex3f(x2, y2, z2);

            rlVertex3f(x1, y1, z1);
            rlVertex3f(x4, y4, z4);
            rlVertex3f(x3, y3, z3);
        }
    }
    rlEnd();
}

static dWorldID world;
static dSpaceID space;
static dJointGroupID contactgroup;

void DebugSpace()
{
    int i, j;

    int n = dSpaceGetNumGeoms(space);

    for (i = 0; i < n; i++)
    {
        dGeomID geom = dSpaceGetGeom(space, i);

        dVector3 pos;
        dQuaternion rotation;

        Matrix mat = MatrixIdentity();

        int geomClass = dGeomGetClass(geom);

        if (dGeomGetClass(geom) != dPlaneClass)
        {
            memcpy(pos, dGeomGetPosition(geom), sizeof(pos));
            dGeomGetQuaternion(geom, rotation);
            getTransformationMatrixFromQuaternion(pos, rotation, &mat);
        }

        rlPushMatrix();

        switch (geomClass)
        {
        case dBoxClass:
        {
            dVector3 sides;
            dGeomBoxGetLengths(geom, sides);
            rlMultMatrixf(MatrixToFloat(mat));

            DrawCube((Vector3){0, 0, 0}, sides[0], sides[1], sides[2], (Color){255, 255, 0, 255});
            DrawCubeWires((Vector3){0, 0, 0}, sides[0], sides[1], sides[2], (Color){255, 0, 0, 255});
            break;
        }

        case dSphereClass:
        {
            dReal radius = dGeomSphereGetRadius(geom);
            rlMultMatrixf(MatrixToFloat(mat));
            DrawSphere((Vector3){0, 0, 0}, radius, (Color){255, 0, 0, 255});
            DrawSphereWires((Vector3){0, 0, 0}, radius, 12, 12, (Color){255, 0, 255, 255});
            break;
        }

        case dCapsuleClass:
        {
            dReal radius, length;
            dGeomCapsuleGetParams(geom, &radius, &length);

            rlMultMatrixf(MatrixToFloat(mat));
            DrawDetailedCapsule(length, radius, (Color){0, 0, 255, 255});
            break;
        }

        case dCylinderClass:
        {
            dReal radius, length;
            dGeomCylinderGetParams(geom, &radius, &length);
            rlMultMatrixf(MatrixToFloat(mat));
            DrawDetailedCylinder(length, radius, (Color){255, 0, 255, 255});
            break;
        }

        case dPlaneClass:
        {
            dVector4 n;
            dMatrix3 R, sides;
            dVector3 pos2;
            dGeomPlaneGetParams(geom, n);
            dRFromZAxis(R, n[0], n[1], n[2]);
            for (j = 0; j < 3; j++)
                pos[j] = n[j] * n[3];

            sides[0] = 200;
            sides[2] = 200;
            sides[1] = 0.001;

            for (j = 0; j < 3; j++)
                pos2[j] = pos[j] + 0.1 * n[j];

            // RenderBox(pos, R, sides);
            DrawCube((Vector3){pos2[0], pos2[1], pos2[2]}, sides[0], sides[1], sides[2], (Color){40, 40, 40, 155});

            break;
        }

        case dRayClass:
        {
            dReal length;
            dVector3 start, dir;
            dGeomRayGet(geom, start, dir);
            length = dGeomRayGetLength(geom);
            Vector3 endPoint = {
                (float)start[0] + dir[0] * length,
                (float)start[1] + dir[1] * length,
                (float)start[2] + dir[2] * length};
            DrawLine3D(
                (Vector3){(float)start[0], (float)start[1], (float)start[2]},
                endPoint,
                RED);
            break;
        }
        }

        rlPopMatrix();
    }
}

#define MAX_CONTACTS 8

static void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
    (void)data;
    int i;

    // exit without doing anything if the two bodies are connected by a joint
    dBodyID b1 = dGeomGetBody(o1);
    dBodyID b2 = dGeomGetBody(o2);
    // if (b1==b2) return;
    if (b1 && b2 && dAreConnectedExcluding(b1, b2, dJointTypeContact))
        return;

    // if (!checkColliding(o1)) return;
    // if (!checkColliding(o2)) return;

    // getting these just so can sometimes be a little bit of a black art!
    dContact contact[MAX_CONTACTS]; // up to MAX_CONTACTS contacts per body-body

    // for (int i = 0; i < MAX_CONTACTS; i++)
    // {
    //     // Modo de contato para permitir rolamento
    //     contact[i].surface.mode = dContactBounce | dContactSoftCFM | dContactApprox1 |
    //                             dContactRolling ;

    //     // Configurar atrito para favorecer rolamento
    //     contact[i].surface.mu = dInfinity;  // Alto atrito para superfícies
    //     contact[i].surface.mu2 = 0;         // Baixo atrito lateral

    //     // Reduzir bounce drasticamente
    //     contact[i].surface.bounce = 0.1;
    //     contact[i].surface.bounce_vel = 0.2;

    //     // Ajustar parâmetros de rolamento
    //     contact[i].surface.rho = 0.1;     // Coeficiente de atrito de rolamento
    //     contact[i].surface.rho2 = 0.1;    // Coeficiente de atrito de rolamento secundário

    //     // Parâmetros para contato mais firme
    //     contact[i].surface.soft_cfm = 0.001;
    //     contact[i].surface.soft_erp = 0.8;

    //     // Adicionar movimento natural
    //     contact[i].surface.motion1 = 0.0;
    //     contact[i].surface.motion2 = 0.0;
    // }

    for (i = 0; i < MAX_CONTACTS; i++)
    {
        contact[i].surface.mode = dContactSlip1 | dContactSlip2 |
                                  dContactSoftERP | dContactSoftCFM | dContactApprox1;
        contact[i].surface.mu = dInfinity; // Alto atrito para superfícies
        contact[i].surface.mu2 = 0;        // Baixo atrito lateral
        contact[i].surface.slip1 = 0.0001;
        contact[i].surface.slip2 = 0.001;
        contact[i].surface.soft_erp = 0.5;
        contact[i].surface.soft_cfm = 0.0003;
        //     parâmetros de rolamento
        contact[i].surface.rho = 0.1;  // Coeficiente de atrito de rolamento
        contact[i].surface.rho2 = 0.1; // Coeficiente de atrito de rolamento secundário
        contact[i].surface.bounce = 0.1;
        contact[i].surface.bounce_vel = 0.2;
    }
    int numc = dCollide(o1, o2, MAX_CONTACTS, &contact[0].geom,
                        sizeof(dContact));
    if (numc)
    {
        dMatrix3 RI;
        dRSetIdentity(RI);
        for (i = 0; i < numc; i++)
        {
            dJointID c =
                dJointCreateContact(world, contactgroup, contact + i);
            dJointAttach(c, b1, b2);
        }
    }
}

struct Vehicle
{
    dBodyID bodies[6];
    dGeomID geoms[6];
    dJointID joints[6];
};

void flipVehicle(Vehicle *car)
{
    const dReal *cp = dBodyGetPosition(car->bodies[0]);
    dBodySetPosition(car->bodies[0], cp[0], cp[1] + 2, cp[2]);

    const dReal *R = dBodyGetRotation(car->bodies[0]);
    dReal newR[16];
    dRFromEulerAngles(newR, 0, -atan2(-R[2], R[0]), 0);
    dBodySetRotation(car->bodies[0], newR);

    dReal wheelOffsets[4][3] = {
        {+1.2, -.6, -1},
        {+1.2, -.6, +1},
        {-1.2, -.6, -1},
        {-1.2, -.6, +1}};

    for (int i = 1; i < 5; i++)
    {
        dVector3 pb;
        dBodyGetRelPointPos(car->bodies[0], wheelOffsets[i - 1][0], wheelOffsets[i - 1][1], wheelOffsets[i - 1][2], pb);
        dBodySetPosition(car->bodies[i], pb[0], pb[1], pb[2]);
    }
}
void updateVehicle(Vehicle *car, float accel, float maxAccelForce, float steer, float steerFactor)
{
    float target;
    target = 0;
    if (fabs(accel) > 0.1)
        target = maxAccelForce;
    // dJointSetHinge2Param( car->joints[0], dParamVel2, -accel );
    // dJointSetHinge2Param( car->joints[1], dParamVel2, accel );
    dJointSetHinge2Param(car->joints[2], dParamVel2, -accel);
    dJointSetHinge2Param(car->joints[3], dParamVel2, accel);

    // dJointSetHinge2Param( car->joints[0], dParamFMax2, target );
    // dJointSetHinge2Param( car->joints[1], dParamFMax2, target );
    dJointSetHinge2Param(car->joints[2], dParamFMax2, target);
    dJointSetHinge2Param(car->joints[3], dParamFMax2, target);

    for (int i = 0; i < 2; i++)
    {
        dReal v = steer - dJointGetHinge2Angle1(car->joints[i]);
        v *= steerFactor;
        dJointSetHinge2Param(car->joints[i], dParamVel, v);
    }

    // float rollForce = 0.5f * (dJointGetHinge2Angle2(car->joints[2]) - dJointGetHinge2Angle2(car->joints[3]));
    // dBodyAddTorque(car->bodies[3], 0, rollForce, 0);
    // dBodyAddTorque(car->bodies[4], 0, -rollForce, 0);
}

float rndf(float min, float max)
{
    return ((float)rand() / (float)(RAND_MAX)) * (max - min) + min;
}

int main(void)
{

    const int screenWidth = 1024;
    const int screenHeight = 720;

    InitWindow(screenWidth, screenHeight, "SMEBotic v0.1 by Luis Santos aka DJOKER");
    SetWindowState(FLAG_WINDOW_RESIZABLE);
    SetTargetFPS(60);
    dInitODE();
    world = dWorldCreate();
    space = dHashSpaceCreate(0);
    contactgroup = dJointGroupCreate(0);

    dWorldSetGravity(world, 0, -9.8, 0);
    dWorldSetAutoDisableFlag(world, 1);
    dWorldSetAutoDisableLinearThreshold(world, 0.05);
    dWorldSetAutoDisableAngularThreshold(world, 0.05);
    dWorldSetAutoDisableSteps(world, 4);

    dGeomID ground = dCreatePlane(space, 0, 1, 0, 0);

    Vector3 carScale = (Vector3){2.5, 0.5, 1.4};
    float wheelRadius = 0.5, wheelWidth = 0.4;

    // car body
    Vehicle car;
    dMass m;
    dMassSetBox(&m, 1, carScale.x, carScale.y, carScale.z); // density
    dMassAdjust(&m, 550);                                   // mass

    car.bodies[0] = dBodyCreate(world);
    dBodySetMass(car.bodies[0], &m);
    dBodySetAutoDisableFlag(car.bodies[0], 0);

    car.geoms[0] = dCreateBox(space, carScale.x, carScale.y, carScale.z);
    dGeomSetBody(car.geoms[0], car.bodies[0]);

    dBodySetPosition(car.bodies[0], 0, 2, 0);

    dGeomID front = dCreateBox(space, 1.5, 0.5, 0.9);
    dGeomSetBody(front, car.bodies[0]);
    dGeomSetOffsetPosition(front, carScale.x / 2 - 1.55, carScale.y / 2 + 0.25, 0);

    car.bodies[5] = dBodyCreate(world);
    dBodySetMass(car.bodies[5], &m);
    dBodySetAutoDisableFlag(car.bodies[5], 0);

    dBodySetPosition(car.bodies[5], 0, 1.5, 0);
    car.geoms[5] = dCreateSphere(space, 0.2);
    dGeomSetBody(car.geoms[5], car.bodies[5]);

    car.joints[5] = dJointCreateFixed(world, 0);
    dJointAttach(car.joints[5], car.bodies[0], car.bodies[5]);
    dJointSetFixed(car.joints[5]);

    // wheels
    dMassSetCylinder(&m, 1, 3, wheelRadius, wheelWidth);
    dMassAdjust(&m, 2); // mass
    dQuaternion q;
    dQFromAxisAndAngle(q, 0, 0, 1, M_PI * 0.5);
    for (int i = 1; i <= 4; ++i)
    {
        car.bodies[i] = dBodyCreate(world);
        dBodySetMass(car.bodies[i], &m);
        dBodySetQuaternion(car.bodies[i], q);
        car.geoms[i] = dCreateCylinder(space, wheelRadius, wheelWidth);
        dGeomSetBody(car.geoms[i], car.bodies[i]);
        dBodySetFiniteRotationMode(car.bodies[i], 1);
        dBodySetAutoDisableFlag(car.bodies[i], 0);
    }

    const dReal *cp = dBodyGetPosition(car.bodies[0]);
    dBodySetPosition(car.bodies[1], cp[0] + 1.2, cp[1] - .5, cp[2] - 1);
    dBodySetPosition(car.bodies[2], cp[0] + 1.2, cp[1] - .5, cp[2] + 1);
    dBodySetPosition(car.bodies[3], cp[0] - 1.2, cp[1] - .5, cp[2] - 1);
    dBodySetPosition(car.bodies[4], cp[0] - 1.2, cp[1] - .5, cp[2] + 1);

    // hinge2 (combined steering / suspension / motor !)
    for (int i = 0; i < 4; ++i)
    {
        car.joints[i] = dJointCreateHinge2(world, 0);
        dJointAttach(car.joints[i], car.bodies[0], car.bodies[i + 1]);
        const dReal *wPos = dBodyGetPosition(car.bodies[i + 1]);
        dJointSetHinge2Anchor(car.joints[i], wPos[0], wPos[1], wPos[2]);

        dReal axis1[] = {0, -1.0, 0};
        dReal axis2[] = {0, 0, ((i % 2) == 0) ? -1.0 : 1.0};

        dJointSetHinge2Axes(car.joints[i], axis1, axis2);
        // dJointSetHinge2Axis1(joints[i], 0, 1, 0);
        // dJointSetHinge2Axis2(joints[i], 0, 0, ((i % 2) == 0) ? -1 : 1);

        dJointSetHinge2Param(car.joints[i], dParamLoStop, 0);
        dJointSetHinge2Param(car.joints[i], dParamHiStop, 0);
        dJointSetHinge2Param(car.joints[i], dParamLoStop, 0);
        dJointSetHinge2Param(car.joints[i], dParamHiStop, 0);
        dJointSetHinge2Param(car.joints[i], dParamFMax, 1500);

        dJointSetHinge2Param(car.joints[i], dParamVel2, dInfinity);
        dJointSetHinge2Param(car.joints[i], dParamFMax2, 1500);

        dJointSetHinge2Param(car.joints[i], dParamSuspensionERP, 0.7);
        dJointSetHinge2Param(car.joints[i], dParamSuspensionCFM, 0.0025);

        // steering
        if (i < 2)
        {
            dJointSetHinge2Param(car.joints[i], dParamFMax, 500);
            dJointSetHinge2Param(car.joints[i], dParamLoStop, -0.5);
            dJointSetHinge2Param(car.joints[i], dParamHiStop, 0.5);
            dJointSetHinge2Param(car.joints[i], dParamLoStop, -0.5);
            dJointSetHinge2Param(car.joints[i], dParamHiStop, 0.5);
            dJointSetHinge2Param(car.joints[i], dParamFudgeFactor, 0.1);
        }
    }
    // disable motor on front wheels
    dJointSetHinge2Param(car.joints[0], dParamFMax2, 0);
    dJointSetHinge2Param(car.joints[1], dParamFMax2, 0);

    Camera camera = {0};
    camera.position = (Vector3){5.0f, 5.0f, 5.0f};
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    bool setCamera = false;

    float accel = 0, steer = 0;
    Vector3 debug = {0};


    bool handbrake = false;
    double frameTime = 0;
    double physTime = 0;
    const double physSlice = 1.0 / 60.0;
    const int maxPsteps = 6;

    while (!WindowShouldClose())
    {
        if (setCamera)
            UpdateCamera(&camera, CAMERA_FREE);

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
        {
            if (setCamera)
                setCamera = false;
            else
                setCamera = true;
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
        {
            setCamera = false;
        }

        if (IsKeyPressed(KEY_ONE))
        {
            dBodyID cubeBody = dBodyCreate(world);
            dMass cubeMass;
            dMassSetBox(&cubeMass, 1.0, 1.0, 1.0, 1.0);
            dBodySetMass(cubeBody, &cubeMass);
            dGeomID cubeGeom = dCreateBox(space, 1.0, 1.0, 1.0);
            dGeomSetBody(cubeGeom, cubeBody);
            dBodySetPosition(cubeBody, 0, 2, 0.0);
        }

        if (IsKeyPressed(KEY_TWO))
        {
            dBodyID sphereBody = dBodyCreate(world);
            dGeomID sphereGeom = dCreateSphere(space, 0.5);
            dMass sphereMass;
            dMassSetSphere(&sphereMass, 1.0, 0.5);
            dBodySetMass(sphereBody, &sphereMass);
            dGeomSetBody(sphereGeom, sphereBody);
            dBodySetPosition(sphereBody, 2, 5, 0);
        }

        if (IsKeyPressed(KEY_THREE))
        {
            dBodyID capsuleBody = dBodyCreate(world);
            dMass capsuleMass;
            dGeomID capsuleGeom = dCreateCapsule(space, 0.3, 1.0);
            dMassSetCapsule(&capsuleMass, 1, 3, 0.3, 1.0);
            dMassAdjust(&capsuleMass, 1.0);
            dBodySetPosition(capsuleBody, -2, 4, 0);
            dGeomSetBody(capsuleGeom, capsuleBody);
        }

        if (IsKeyPressed(KEY_FOUR))
        {
            dBodyID cylinderBody = dBodyCreate(world);
            dReal r = 1.0;
            dReal l = 0.4;
            dGeomID cylinderGeom = dCreateCylinder(space, r, l);
            dMass cylinderMass;
            int dir = 1; // (1=x, 2=y, 3=z)
            dMassSetCylinder(&cylinderMass, 1, dir, r, l);
            dBodySetMass(cylinderBody, &cylinderMass);
            dMassAdjust(&cylinderMass, 1.0);
            dBodySetPosition(cylinderBody, -2, 4, 0);
            dGeomSetBody(cylinderGeom, cylinderBody);

            // // Configurar bits de colisão para permitir colisão entre cilindros
            // unsigned long category = 0x0001;    // Categoria para cilindros
            // unsigned long collide = 0xFFFFFFFF; // Colide com todas as categorias
            // dGeomSetCategoryBits(cylinderGeom, category);
            // dGeomSetCollideBits(cylinderGeom, collide);
        }

        if (IsKeyDown(KEY_SPACE))
        {
            dJointSetHinge2Param(car.joints[2], dParamVel2, 0);
            dJointSetHinge2Param(car.joints[3], dParamVel2, 0);

            dJointSetHinge2Param(car.joints[2], dParamFMax2, 10000); 
            dJointSetHinge2Param(car.joints[3], dParamFMax2, 10000);
        }
        else
        {
         
            dJointSetHinge2Param(car.joints[2], dParamFMax2, 1500);
            dJointSetHinge2Param(car.joints[3], dParamFMax2, 1500);
        }

        if (IsKeyPressed(KEY_R))
            flipVehicle(&car);

        accel *= .99;
        if (IsKeyDown(KEY_UP))
            accel += 2;
        if (IsKeyDown(KEY_DOWN))
            accel -= 2;
        if (accel > 50)
            accel = 50;
        if (accel < -15)
            accel = -15;

        if (IsKeyDown(KEY_RIGHT))
            steer -= .1;
        if (IsKeyDown(KEY_LEFT))
            steer += .1;
        if (!IsKeyDown(KEY_RIGHT) && !IsKeyDown(KEY_LEFT))
            steer *= .5;
        if (steer > .5)
            steer = .5;
        if (steer < -.5)
            steer = -.5;

        updateVehicle(&car, accel, 800.0, steer, 10.0);

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        DrawGrid(10, 1.0f);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){1.0f, 0.5f, 0.0f}, RED);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){0.0f, 0.5f, 1.0f}, GREEN);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){0.0f, 1.5f, 0.0f}, BLUE);

        frameTime += GetFrameTime();
        int pSteps = 0;
        physTime = GetTime();

        while (frameTime > physSlice)
        {

            dSpaceCollide(space, 0, &nearCallback);

            // step the world
            dWorldQuickStep(world, physSlice);
            dJointGroupEmpty(contactgroup);

            frameTime -= physSlice;
            pSteps++;
            if (pSteps > maxPsteps)
            {
                frameTime = 0;
                break;
            }
        }

        physTime = GetTime() - physTime;

        DebugSpace();

        EndMode3D();

        if (pSteps > maxPsteps)
            DrawText("WARNING CPU overloaded lagging real time", 10, 0, 20, RED);
        DrawText(TextFormat("%2i FPS", GetFPS()), 10, 20, 20, WHITE);
        DrawText(TextFormat("accel %4.4f", accel), 10, 40, 20, WHITE);
        DrawText(TextFormat("steer %4.4f", steer), 10, 60, 20, WHITE);
         DrawText(TextFormat("debug %4.4f %4.4f %4.4f", debug.x, debug.y, debug.z), 10, 100, 20, WHITE);
        DrawText(TextFormat("Phys steps per frame %i", pSteps), 10, 120, 20, WHITE);
        DrawText(TextFormat("Phys time per frame %i", physTime), 10, 140, 20, WHITE);
        DrawText(TextFormat("total time per frame %i", frameTime), 10, 160, 20, WHITE);

        //  DrawText(TextFormat("roll %.4f", fabs(roll)), 10, 200, 20, WHITE);

        const double *v = dBodyGetLinearVel(car.bodies[0]);
        float vel = Vector3Length((Vector3){v[0], v[1], v[2]}) * 2.23693629f;
        DrawText(TextFormat("mph %.4f", vel), 10, 220, 20, WHITE);

        EndDrawing();
    }

    dJointGroupDestroy(contactgroup);
    dGeomDestroy(ground);
    dSpaceDestroy(space);
    dWorldDestroy(world);
    dCloseODE();

    CloseWindow();

    return 0;
}
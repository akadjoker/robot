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

        // Face lateral - segundo triângulo
        rlVertex3f(x2, y2, -halfLength); 
        rlVertex3f(x2, y2, halfLength); 
        rlVertex3f(x1, y1, halfLength);  

        // Face superior
        rlColor4ub(color.r * 0.9, color.g * 0.9, color.b * 0.9, color.a);
        rlVertex3f(0, 0, halfLength);   
        rlVertex3f(x1, y1, halfLength); 
        rlVertex3f(x2, y2, halfLength); 

        // Face inferior
        rlColor4ub(color.r * 0.8, color.g * 0.8, color.b * 0.8, color.a);
        rlVertex3f(0, 0, -halfLength);   
        rlVertex3f(x2, y2, -halfLength); 
        rlVertex3f(x1, y1, -halfLength); 
    }
    rlEnd();
}

void DrawDetailedCapsule(float length, float radius, Color color)
{
    const int segments = 24;
    const int rings = 12;    //  semiesferas
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
           // DrawCylinderWires((Vector3){0, 0, 0}, radius, radius, length, 20, (Color){255, 255, 0, 255});
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

#define MAX_CONTACTS 20



void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
    dBodyID b1 = dGeomGetBody(o1);
    dBodyID b2 = dGeomGetBody(o2);
    
    if (b1 && b2 && dAreConnectedExcluding(b1, b2, dJointTypeContact))
        return;

    dContact contact[MAX_CONTACTS];
    for (int i = 0; i < MAX_CONTACTS; i++)
    {
        // Modo de contato para permitir rolamento
        contact[i].surface.mode = dContactBounce | dContactSoftCFM | dContactApprox1 | 
                                dContactRolling ;
        
        // Configurar atrito para favorecer rolamento
        contact[i].surface.mu = dInfinity;  // Alto atrito para superfícies
        contact[i].surface.mu2 = 0;         // Baixo atrito lateral
        
        // Reduzir bounce drasticamente
        contact[i].surface.bounce = 0.1;
        contact[i].surface.bounce_vel = 0.2;
        
        // Ajustar parâmetros de rolamento
        contact[i].surface.rho = 0.1;     // Coeficiente de atrito de rolamento
        contact[i].surface.rho2 = 0.1;    // Coeficiente de atrito de rolamento secundário
        
        // Parâmetros para contato mais firme
        contact[i].surface.soft_cfm = 0.001;
        contact[i].surface.soft_erp = 0.8;
        
        // Adicionar movimento natural
        contact[i].surface.motion1 = 0.0;
        contact[i].surface.motion2 = 0.0;
    }

    int numc = dCollide(o1, o2, MAX_CONTACTS, &contact[0].geom, sizeof(dContact));
    if (numc > 0)
    {
        for (int i = 0; i < numc; i++)
        {
            dJointID c = dJointCreateContact(world, contactgroup, &contact[i]);
            dJointAttach(c, b1, b2);
            
            // Aplicar força de rolamento se for uma esfera
            // if (dGeomGetClass(o1) == dSphereClass || dGeomGetClass(o2) == dSphereClass)
            // {
            //     dBodyID sphereBody = (dGeomGetClass(o1) == dSphereClass) ? b1 : b2;
            //     const dReal* rot = dBodyGetAngularVel(sphereBody);
            //     //dBodyAddTorque(sphereBody, -rot[0] * 0.1, -rot[1] * 0.1, -rot[2] * 0.1);
            //     dBodyAddTorque(sphereBody, 0,0,1);
            // }
        }
    }
}


int main(void)
{

    const int screenWidth = 1024;
    const int screenHeight = 720;

    InitWindow(screenWidth, screenHeight, "SMEBotic v0.1 by Luis Santos aka DJOKER");
    SetTargetFPS(60);
    dInitODE();
    world = dWorldCreate();
    space = dHashSpaceCreate(0);
    contactgroup = dJointGroupCreate(0);

    dWorldSetGravity(world, 0, -9.8, 0);

    // dWorldSetCFM(world, 1e-5);

    // dWorldSetLinearDamping(world, 0.00001);
    // dWorldSetAngularDamping(world, 0.005);
    // dWorldSetMaxAngularSpeed(world, 200);

    // dWorldSetContactSurfaceLayer(world, 0.001);

    dGeomID ground = dCreatePlane(space, 0, 1, 0, 0);

    dBodyID cubeBody = dBodyCreate(world);
    dMass cubeMass;
    dMassSetBox(&cubeMass, 1.0, 1.0, 1.0, 1.0);
    dMassAdjust(&cubeMass, 1.0);
    dBodySetMass(cubeBody, &cubeMass);
    dGeomID cubeGeom = dCreateBox(space, 1.0, 1.0, 1.0);
    dGeomSetBody(cubeGeom, cubeBody);
    dBodySetPosition(cubeBody, -2, 5, 0);

    dBodyID sphereBody = dBodyCreate(world);
    dMass sphereMass;
    dMassSetSphere(&sphereMass, 1.0, 0.5);
    // dMassAdjust(&sphereMass, 1.0);
    dBodySetMass(sphereBody, &sphereMass);
    dGeomID sphereGeom = dCreateSphere(space, 0.5);
    dGeomSetBody(sphereGeom, sphereBody);
    dBodySetPosition(sphereBody, 2, 5, 0);

    dBodyID capsuleBody = dBodyCreate(world);
    dMass capsuleMass;
    dGeomID capsuleGeom = dCreateCapsule(space, 0.3, 1.0);
    dMassSetCapsule(&capsuleMass, 1, 3, 0.3, 1.0);
    // dMassAdjust(&capsuleMass, 1.0);
    dBodySetPosition(capsuleBody, -2, 4, 0);
    dGeomSetBody(capsuleGeom, capsuleBody);

    dBodyID cylinderBody = dBodyCreate(world);
    dReal r = 0.9;
    dReal l = 0.6;
    dGeomID cylinderGeom = dCreateCylinder(space, r, l);
    dMass cylinderMass;
    // (1=x, 2=y, 3=z).
    int dir = 1;
    dMassSetCylinder(&cylinderMass, 1, dir, r, l);
    dBodySetMass(cylinderBody, &cylinderMass);

    dBodySetPosition(cylinderBody, 0, 6, 2);
    dGeomSetBody(cylinderGeom, cylinderBody);
    // dMassAdjust(&cylinderMass, 1.0);
    // dMassSetCylinderTotal(&cylinderMass, 1, 3, 0.3, 1.2);

    // dGeomSetCategoryBits(ground, 1ul);
    // dGeomSetCategoryBits(cubeGeom, 2ul);

    // dGeomSetCollideBits(ground, ~2ul);
    // dGeomSetCollideBits(cubeGeom, ~1ul);

    Camera camera = {0};
    camera.position = (Vector3){10.0f, 10.0f, 10.0f};
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    bool setCamera = false;

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
            // dMassAdjust(&cubeMass, 1.0);
            dBodySetMass(cubeBody, &cubeMass);
            dGeomID cubeGeom = dCreateBox(space, 1.0, 1.0, 1.0);
            dGeomSetBody(cubeGeom, cubeBody);
            dBodySetPosition(cubeBody, -2, 5, 0);
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
            dBodyAddTorque(sphereBody, 1,0,1);
            //dBodySetAngularDamping(sphereBody, 0.01);
            //sBodySetLinearDamping(sphereBody, 0.02);
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

            // Configurar bits de colisão para permitir colisão entre cilindros
            unsigned long category = 0x0001;    // Categoria para cilindros
            unsigned long collide = 0xFFFFFFFF; // Colide com todas as categorias
            dGeomSetCategoryBits(cylinderGeom, category);
            dGeomSetCollideBits(cylinderGeom, collide);
        }

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        DrawGrid(10, 1.0f);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){1.0f, 0.5f, 0.0f}, RED);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){0.0f, 0.5f, 1.0f}, GREEN);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){0.0f, 1.5f, 0.0f}, BLUE);
        dSpaceCollide(space, 0, &nearCallback);
        ///        dWorldStep(world, 0.05);
        const dReal stepsize = 1.0 / 60.0;
        dWorldQuickStep(world, stepsize);
        dJointGroupEmpty(contactgroup);

        DebugSpace();

        EndMode3D();

        EndDrawing();
    }

    dGeomDestroy(cubeGeom);
    dGeomDestroy(sphereGeom);
    dBodyDestroy(cubeBody);
    dBodyDestroy(sphereBody);
    dJointGroupDestroy(contactgroup);
    dGeomDestroy(ground);
    dSpaceDestroy(space);
    dWorldDestroy(world);
    dCloseODE();

    CloseWindow();

    return 0;
}
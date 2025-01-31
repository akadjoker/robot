#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
#include "ode/ode.h"
#include "RigidBody.hpp"
#include "Joint.hpp"
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/attr.h>
#include <sys/stat.h>
#include <cstdint>
#include <string>

namespace py = pybind11;

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

dWorldID world;
dSpaceID space;
dJointGroupID contactgroup;

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

PYBIND11_EMBEDDED_MODULE(sme, m)
{
    m.doc() = "Ray Lib with ODE";

    py::class_<Vector2>(m, "Vec2")
        .def(py::init<>())             // Construtor vazio
        .def(py::init<float, float>()) // Construtor com valores x e y
        .def_readwrite("x", &Vector2::x)
        .def_readwrite("y", &Vector2::y);

    py::class_<Vector3>(m, "Vec3")
        .def(py::init<>())                    // Construtor vazio
        .def(py::init<float, float, float>()) // Construtor com valores x, y e z
        .def_readwrite("x", &Vector3::x)
        .def_readwrite("y", &Vector3::y)
        .def_readwrite("z", &Vector3::z);

    py::class_<RigidBody>(m, "RigidBody")
        .def(py::init<bool>(), py::arg("isStatic") = false)
        .def("SetPosition", &RigidBody::setPosition, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("SetBodyPosition", &RigidBody::setBodyPosition, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("SetGeomPosition", &RigidBody::setGeomPosition, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("SetId", &RigidBody::setId, py::arg("id"))
        .def("GetId", &RigidBody::getId);

    py::class_<Cube, RigidBody, std::unique_ptr<Cube, py::nodelete>>(m, "Cube")
        .def(py::init<dReal, dReal, dReal, dReal, bool>(),
             py::arg("width"), py::arg("height"), py::arg("depth"),
             py::arg("mass"), py::arg("isStatic") = false);
    py::class_<Sphere, RigidBody, std::unique_ptr<Sphere, py::nodelete>>(m, "Sphere")
        .def(py::init<dReal, dReal, bool>(),
             py::arg("radius"), py::arg("mass"), py::arg("isStatic") = false);

    py::class_<Capsule, RigidBody, std::unique_ptr<Capsule, py::nodelete>>(m, "Capsule")
        .def(py::init<int, dReal, dReal, dReal, bool>(),
             py::arg("direction"), py::arg("radius"), py::arg("length"),
             py::arg("mass"), py::arg("isStatic") = false);

    py::class_<Cylinder, RigidBody, std::unique_ptr<Cylinder, py::nodelete>>(m, "Cylinder")
        .def(py::init<int, dReal, dReal, dReal, bool>(),
             py::arg("direction"), py::arg("radius"), py::arg("length"),
             py::arg("mass"), py::arg("isStatic") = false);

    // Vehicle

    py::class_<Vehicle, RigidBody, std::unique_ptr<Vehicle, py::nodelete>>(m, "Vehicle")
        .def(py::init<dReal, dReal, dReal, dReal>(), py::arg("mass"), py::arg("width"), py::arg("height"), py::arg("depth"))
        .def("AddWheel", &Vehicle::addWheel, py::arg("radius"), py::arg("width"), py::arg("mass"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("ax"), py::arg("ay"), py::arg("az"))
        .def("SetSuspension", &Vehicle::setSuspension, py::arg("index"), py::arg("ERP"), py::arg("CFM"))

        .def("SetParam", &Vehicle::SetParam, py::arg("index"), py::arg("type"), py::arg("valeu"))
        .def("SetAcceleration", &Vehicle::setAcceleration, py::arg("index"), py::arg("accel"), py::arg("maxAccelForce"))
        .def("Break", &Vehicle::applyBrake, py::arg("index"), py::arg("force"))
        .def("SetSteering", &Vehicle::setSteering, py::arg("index"), py::arg("steer"), py::arg("steerFactor"));

    py::class_<Joint>(m, "Joint")
        .def(py::init<>())
        .def("SetId", &Joint::setId, py::arg("id"))
        .def("GetId", &Joint::getId)
        .def("Attach", &Joint::attach, py::arg("body1"), py::arg("body2"));

    py::class_<HingeJoint, Joint, std::unique_ptr<HingeJoint, py::nodelete>>(m, "HingeJoint")
        .def(py::init<>())
        .def("setAxis", &HingeJoint::setAxis, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("setAnchor", &HingeJoint::setAnchor, py::arg("x"), py::arg("y"), py::arg("z"));

    py::class_<Hinge2Joint, Joint, std::unique_ptr<Hinge2Joint, py::nodelete>>(m, "Hinge2Joint")
        .def(py::init<>())
        .def("setAnchor", &Hinge2Joint::setAnchor, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("setAxes", &Hinge2Joint::setAxes, py::arg("ax"), py::arg("ay"), py::arg("az"), py::arg("bx"), py::arg("by"), py::arg("bz"))
        .def("setSuspension", &Hinge2Joint::setSuspension, py::arg("ERP"), py::arg("CFM"));

    py::class_<FixedJoint, Joint, std::unique_ptr<FixedJoint, py::nodelete>>(m, "FixedJoint")
        .def(py::init<>())
        .def("setFixed", &FixedJoint::setFixed);

    py::class_<SliderJoint, Joint, std::unique_ptr<SliderJoint, py::nodelete>>(m, "SliderJoint")
        .def(py::init<>())
        .def("setAxis", &SliderJoint::setAxis, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("setLimits", &SliderJoint::setLimits, py::arg("lo"), py::arg("hi"));

    py::class_<BallJoint, Joint, std::unique_ptr<BallJoint, py::nodelete>>(m, "BallJoint")
        .def(py::init<>())
        .def("setAnchor", &BallJoint::setAnchor, py::arg("x"), py::arg("y"), py::arg("z"));

    py::class_<AMotorJoint, Joint, std::unique_ptr<AMotorJoint, py::nodelete>>(m, "MotorJoint")
        .def(py::init<>())
        .def("setNumAxes", &AMotorJoint::setNumAxes, py::arg("num"))
        .def("setAxis", &AMotorJoint::setAxis, py::arg("axisIndex"), py::arg("relativeTo"), py::arg("x"), py::arg("y"), py::arg("z"))
        .def("setVelocity", &AMotorJoint::setVelocity, py::arg("axisIndex"), py::arg("velocity"))
        .def("setFMax", &AMotorJoint::setFMax, py::arg("axisIndex"), py::arg("maxForce"));

    // world

    m.def("ClearWorld", &ClearWorld);
    m.def("SetGravity", [](dReal x, dReal y, dReal z)
          { dWorldSetGravity(world, x, y, z); });

    m.def("SetAutoDisableFlag", [](bool flag)
          { dWorldSetAutoDisableFlag(world, flag); });

    m.def("SetAutoDisableLinearThreshold", [](dReal threshold)
          { dWorldSetAutoDisableLinearThreshold(world, threshold); });

    m.def("SetAutoDisableAngularThreshold", [](dReal threshold)
          { dWorldSetAutoDisableAngularThreshold(world, threshold); });

    m.def("SetAutoDisableSteps", [](int steps)
          { dWorldSetAutoDisableSteps(world, steps); });

    // keys
    m.def("IsKeyPressed", &IsKeyPressed, py::arg("key"));
    m.def("IsKeyReleased", &IsKeyReleased, py::arg("key"));
    m.def("IsKeyDown", &IsKeyDown, py::arg("key"));
    m.def("IsKeyUp", &IsKeyUp, py::arg("key"));

    // mouse
    m.def("IsMouseButtonPressed", &IsMouseButtonPressed, py::arg("button"));
    m.def("IsMouseButtonReleased", &IsMouseButtonReleased, py::arg("button"));
    m.def("IsMouseButtonDown", &IsMouseButtonDown, py::arg("button"));
    m.def("IsMouseButtonUp", &IsMouseButtonUp, py::arg("button"));
    m.def("GetMousePosition", []
          { return GetMousePosition(); });
    m.def("GetMouseDelta", &GetMouseDelta);

    // draw
    m.def("DrawFPS", &DrawFPS, py::arg("posX"), py::arg("posY"));
}

#define MAX_CONTACTS 8

// Enum para tipos de superfície
enum SurfaceType
{
    ASPHALT,   // Asfalto
    CONCRETE,  // Concreto
    DIRT_ROAD, // Estrada de terra
    GRASS,     // Grama
    ICE,       // Gelo
    SAND       // Areia
};

// Estrutura para configurações de contato do veículo
struct VehicleContactParams
{
    float mu;         // Coeficiente de atrito principal
    float mu2;        // Coeficiente de atrito secundário
    float slip1;      // Deslizamento direção 1
    float slip2;      // Deslizamento direção 2
    float soft_erp;   // Erro de redução
    float soft_cfm;   // Mistura de força de restrição
    float rho;        // Atrito de rolamento
    float bounce;     // Elasticidade
    float bounce_vel; // Velocidade mínima para ressalto
};

// Função para obter parâmetros de contato baseado no tipo de superfície
VehicleContactParams getVehicleContactParams(SurfaceType surface, bool isHighPerformance = false)
{
    VehicleContactParams params;

    switch (surface)
    {
    case ASPHALT:
        params.mu = isHighPerformance ? 1.2f : 0.8f; // Atrito alto para corrida
        params.mu2 = 0.3f;
        params.slip1 = 0.0001f;
        params.slip2 = 0.001f;
        params.soft_erp = 0.6f;
        params.soft_cfm = 0.0002f;
        params.rho = 0.05f;
        params.bounce = 0.1f;
        params.bounce_vel = 0.2f;
        break;

    case CONCRETE:
        params.mu = 0.9f;
        params.mu2 = 0.2f;
        params.slip1 = 0.0002f;
        params.slip2 = 0.002f;
        params.soft_erp = 0.5f;
        params.soft_cfm = 0.0003f;
        params.rho = 0.1f;
        params.bounce = 0.05f;
        params.bounce_vel = 0.1f;
        break;

    case DIRT_ROAD:
        params.mu = 0.6f;
        params.mu2 = 0.4f;
        params.slip1 = 0.001f;
        params.slip2 = 0.005f;
        params.soft_erp = 0.4f;
        params.soft_cfm = 0.0005f;
        params.rho = 0.2f;
        params.bounce = 0.2f;
        params.bounce_vel = 0.3f;
        break;

    case GRASS:
        params.mu = 0.4f;
        params.mu2 = 0.5f;
        params.slip1 = 0.002f;
        params.slip2 = 0.01f;
        params.soft_erp = 0.3f;
        params.soft_cfm = 0.001f;
        params.rho = 0.3f;
        params.bounce = 0.3f;
        params.bounce_vel = 0.4f;
        break;

    case ICE:
        params.mu = 0.1f;
        params.mu2 = 0.05f;
        params.slip1 = 0.005f;
        params.slip2 = 0.02f;
        params.soft_erp = 0.2f;
        params.soft_cfm = 0.002f;
        params.rho = 0.05f;
        params.bounce = 0.4f;
        params.bounce_vel = 0.5f;
        break;

    case SAND:
        params.mu = 0.5f;
        params.mu2 = 0.6f;
        params.slip1 = 0.003f;
        params.slip2 = 0.015f;
        params.soft_erp = 0.3f;
        params.soft_cfm = 0.0008f;
        params.rho = 0.4f;
        params.bounce = 0.2f;
        params.bounce_vel = 0.3f;
        break;

    default:
        // Configuração padrão (asfalto)
        params.mu = 0.8f;
        params.mu2 = 0.3f;
        params.slip1 = 0.0001f;
        params.slip2 = 0.001f;
        params.soft_erp = 0.5f;
        params.soft_cfm = 0.0003f;
        params.rho = 0.1f;
        params.bounce = 0.1f;
        params.bounce_vel = 0.2f;
        break;
    }

    return params;
}

// Função para aplicar parâmetros de contato
void applyVehicleContactParams(dContact &contact, const VehicleContactParams &params)
{
    contact.surface.mode = dContactSlip1 | dContactSlip2 |
                           dContactSoftERP | dContactSoftCFM |
                           dContactApprox1;

    contact.surface.mu = params.mu;
    contact.surface.mu2 = params.mu2;
    contact.surface.slip1 = params.slip1;
    contact.surface.slip2 = params.slip2;
    contact.surface.soft_erp = params.soft_erp;
    contact.surface.soft_cfm = params.soft_cfm;
    contact.surface.rho = params.rho;
    contact.surface.bounce = params.bounce;
    contact.surface.bounce_vel = params.bounce_vel;
}

static void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
    (void)data;
    int i;

    // exit without doing anything if the two bodies are connected by a joint
    dBodyID b1 = dGeomGetBody(o1);
    dBodyID b2 = dGeomGetBody(o2);
    if (b1 == b2)
        return;
    if (b1 && b2 && dAreConnectedExcluding(b1, b2, dJointTypeContact))
        return;

    // if (!checkColliding(o1)) return;
    // if (!checkColliding(o2)) return;

    // getting these just so can sometimes be a little bit of a black art!
    dContact contact[MAX_CONTACTS];

    for (i = 0; i < MAX_CONTACTS; i++)
    {
          VehicleContactParams params = getVehicleContactParams(ASPHALT, true);
         applyVehicleContactParams(contact[i], params);
       
        // contact[i].surface.mode = dContactSlip1 | dContactSlip2 |
        //                           dContactSoftERP | dContactSoftCFM | dContactApprox1;
        // contact[i].surface.mu = 0.8; // Alto atrito para superfícies
        // contact[i].surface.mu2 = 0;  // Baixo atrito lateral
        // contact[i].surface.slip1 = 0.0001;
        // contact[i].surface.slip2 = 0.001;
        // contact[i].surface.soft_erp = 0.5;
        // contact[i].surface.soft_cfm = 0.0003;
        // //     parâmetros de rolamento
        // contact[i].surface.rho = 0.1;  // Coeficiente de atrito de rolamento
        // contact[i].surface.rho2 = 0.1; // Coeficiente de atrito de rolamento secundário
        // contact[i].surface.bounce = 0.1;
        // contact[i].surface.bounce_vel = 0.2;
    }
    int numc = dCollide(o1, o2, MAX_CONTACTS, &contact[0].geom, sizeof(dContact));
    if (numc)
    {
        for (i = 0; i < numc; i++)
        {
            dJointID c = dJointCreateContact(world, contactgroup, contact + i);
            dJointAttach(c, b1, b2);
        }
    }
}

#define DEF_LOAD 0
#define DEF_UNLOAD 1
#define DEF_LOOP 2

bool HotRealod = false;
int WIDTH, HEIGHT;
bool FULLSCREEN = false;
int X, Y;
std::string TITLE;

py::object load_function;
py::object unload_function;
py::object script_loop_function;

bool Exists[3] = {false, false, false};

uint64_t GetFileModifiedTime(const std::string &filePath)
{
    struct stat fileStat;
    if (stat(filePath.c_str(), &fileStat) == 0)
    {
        return static_cast<uint64_t>(fileStat.st_mtime);
    }
    return 0; // Retorna 0 se o arquivo não existir ou ocorrer um erro
}
bool load_script(const std::string &script)
{
    try
    {
        py::exec("if 'HOT_RELOAD' not in globals(): HOT_RELOAD = False");
        py::exec("if 'WINDOW_WIDTH' not in globals(): WINDOW_WIDTH = 800");
        py::exec("if 'WINDOW_HEIGHT' not in globals(): WINDOW_HEIGHT = 600");
        py::exec("if 'WINDOW_FULLSCREEN' not in globals(): WINDOW_FULLSCREEN = False");
        py::exec("if 'WINDOW_TITLE' not in globals(): WINDOW_TITLE = 'SME Robotics Engine by Luis Santos AKA DJOKER'");
        py::exec("if 'WINDOW_X' not in globals(): WINDOW_X = 0");
        py::exec("if 'WINDOW_Y' not in globals(): WINDOW_Y = 0");

        py::exec("global load, unload, loop");

        py::exec("load = globals().get('load')");
        py::exec("unload = globals().get('unload')");
        py::exec("loop = globals().get('loop')");

        py::eval_file(script);

        WIDTH = py::globals()["WINDOW_WIDTH"].cast<int>();
        HEIGHT = py::globals()["WINDOW_HEIGHT"].cast<int>();
        FULLSCREEN = py::globals()["WINDOW_FULLSCREEN"].cast<bool>();
        TITLE = py::globals()["WINDOW_TITLE"].cast<std::string>();
        X = py::globals()["WINDOW_X"].cast<int>();
        Y = py::globals()["WINDOW_Y"].cast<int>();

        // Cache functions
        load_function = py::globals()["load"];
        unload_function = py::globals()["unload"];
        script_loop_function = py::globals()["loop"];

        HotRealod = py::globals()["HOT_RELOAD"].cast<bool>();

        Exists[DEF_LOAD] = !load_function.is_none();
        Exists[DEF_UNLOAD] = !unload_function.is_none();
        Exists[DEF_LOOP] = !script_loop_function.is_none();
    }
    catch (py::error_already_set &e)
    {
        // Erros Python
        PyErr_Print();
        return false;
    }
    catch (const std::exception &e)
    {
        // Erros padrão do C++
        std::string err = "Uncaught std::exception: ";
        err += e.what();
        printf("%s\n", err.c_str());
        return false;
    }
    catch (...)
    {
        // Erros genéricos
        std::string err = "Uncaught exception occurred!";
        printf("%s\n", err.c_str());
        return false;
    }

    return true;
}

void load_def()
{
    try
    {
        if (Exists[DEF_LOAD])
        {
            load_function();
        }
    }
    catch (py::error_already_set &e)
    {
        PyErr_Print();
    }
}

void unload_script()
{
    try
    {
        if (Exists[DEF_UNLOAD])
        {
            unload_function();
        }
        // py::exec("import gc; gc.collect()");
        // py::exec("globals().clear()");
    }
    catch (py::error_already_set &e)
    {
        PyErr_Print();
    }
}

void script_loop()
{
    try
    {
        if (Exists[DEF_LOOP])
        {
            script_loop_function();
        }
    }
    catch (py::error_already_set &e)
    {
        PyErr_Print();
    }
}

int main(int argc, char *argv[])
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

    std::string script = "main.py";
    if (argc >= 2)
    {
        script = argv[1];
    }
    if (!FileExists(script.c_str()))
    {
        printf("Script '%s' not found!\n", script.c_str());
        dJointGroupDestroy(contactgroup);
        dSpaceDestroy(space);
        dWorldDestroy(world);
        dCloseODE();
        CloseWindow();
        return 1;
    }
    uint64_t last_modified_time = GetFileModTime(script.c_str());
    py::scoped_interpreter guard{};
    if (!load_script(script))
    {
        printf("Failed to load script '%s'\n", script.c_str());
        dJointGroupDestroy(contactgroup);
        dSpaceDestroy(space);
        dWorldDestroy(world);
        dCloseODE();
        CloseWindow();
        return 1;
    }
    printf("Script '%s' loaded\n", script.c_str());

    dWorldSetGravity(world, 0, -9.8, 0);
    dWorldSetAutoDisableFlag(world, 1);
    dWorldSetAutoDisableLinearThreshold(world, 0.05);
    dWorldSetAutoDisableAngularThreshold(world, 0.05);
    dWorldSetAutoDisableSteps(world, 4);

    dGeomID ground = dCreatePlane(space, 0, 1, 0, 0);

    double frameTime = 0;
    double physTime = 0;
    const double physSlice = 1.0 / 60.0;
    const int maxPsteps = 6;

    Camera camera = {0};
    camera.position = (Vector3){5.0f, 5.0f, 5.0f};
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    bool setCamera = false;

    load_def();

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

        // if (IsKeyPressed(KEY_ONE))
        // {
        //     dBodyID cubeBody = dBodyCreate(world);
        //     dMass cubeMass;
        //     dMassSetBox(&cubeMass, 1.0, 1.0, 1.0, 1.0);
        //     dBodySetMass(cubeBody, &cubeMass);
        //     dGeomID cubeGeom = dCreateBox(space, 1.0, 1.0, 1.0);
        //     dGeomSetBody(cubeGeom, cubeBody);
        //     dBodySetPosition(cubeBody, 0, 2, 0.0);
        // }

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
        }

        if (HotRealod)
        {

            uint64_t current_modified_time = GetFileModTime(script.c_str());
            if (current_modified_time > last_modified_time)
            {
                printf("Reloading script...\n");
                unload_script();
                if (load_script(script))
                {
                    printf("Script reloaded successfully.\n");
                    load_def();
                    last_modified_time = current_modified_time;
                }
            }
        }

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

        script_loop();
        // DrawFPS(10, 10);

        EndDrawing();
    }

    unload_script();
    ClearWorld();
    dJointGroupDestroy(contactgroup);
    dGeomDestroy(ground);
    dSpaceDestroy(space);
    dWorldDestroy(world);
    dCloseODE();
    CloseWindow();

    return 0;
}
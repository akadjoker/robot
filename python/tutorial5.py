#!/usr/bin/env python

import sys, random, time
from math import *
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import ode

# 🏎️ Variáveis de controle do carro
speed = 0.0  # Velocidade
steer = 0.0  # Direção

# ⚙️ Inicializar o ODE
world = ode.World()
world.setGravity((0, -9.81, 0))
world.setERP(0.8)
world.setCFM(1E-5)
space = ode.Space()
contactgroup = ode.JointGroup()

# Criar um plano para impedir que os objetos caiam infinitamente
floor = ode.GeomPlane(space, (0, 1, 0), 0)

# Listas para armazenar corpos e juntas
bodies = []
geoms = []
wheels = []
joints = []

# 📌 Criar o chassis do carro
def create_chassis(world, space):
    body = ode.Body(world)
    M = ode.Mass()
    M.setBox(500, 1.5, 0.5, 3)  # Massa e tamanho do chassis
    body.setMass(M)
    body.shape = "box"
    body.boxsize = (1.5, 0.5, 3)
    geom = ode.GeomBox(space, lengths=body.boxsize)
    geom.setBody(body)
    return body, geom

# 🛞 Criar uma roda
def create_wheel(world, space):
    body = ode.Body(world)
    M = ode.Mass()
    M.setCylinder(50, 3, 0.3, 0.15)
    body.setMass(M)
    body.shape = "cylinder"
    body.cyl_params = (0.3, 0.15)
    geom = ode.GeomCylinder(space, 0.3, 0.15)
    geom.setBody(body)
    return body, geom

def make_orthogonal(axis1, axis2):
    """Garante que axis2 é ortogonal a axis1"""
    if abs(sum(a * b for a, b in zip(axis1, axis2))) > 1e-5:  # Se não for ortogonal
        return (axis2[2], axis2[0], axis2[1])  # Troca os valores para garantir ortogonalidade
    return axis2

def attach_wheel(world, car_body, wheel_body, x, y, z, is_front):
    """Cria uma junta Hinge2 para conectar uma roda ao chassis"""
    joint = ode.Hinge2Joint(world)
    joint.attach(car_body, wheel_body)
    joint.setAnchor((x, y, z))

  
    axis1 = (0, 1, 0)
    axis2 = (1, 0, 0) 
    #if is_front else (0, 0, 1)
    axis2 = make_orthogonal(axis1, axis2)  # Garante que são perpendiculares

    print(f"🚗 Rodando ({x}, {y}, {z}) -> Axis1: {axis1}, Axis2: {axis2}")
    joint.setParam(ode.ParamSuspensionERP, 0.4)
    joint.setParam(ode.ParamSuspensionCFM, 0.8)

    # 🚗 Travar as rodas traseiras para que não virem
    if not is_front:
        joint.setParam(ode.ParamLoStop, 0)
        joint.setParam(ode.ParamHiStop, 0)

    return joint




# 🚗 Criar o carro completo
def create_car(world, space):
    global wheels, joints
    car_body, car_geom = create_chassis(world, space)
    wheel_positions = [
    (0.8, 0, 1.4, True),  # Frente direita
    (-0.8, 0, 1.4, True), # Frente esquerda
    (0.8, 0, -1.4, False), # Trás direita
    (-0.8, 0, -1.4, False) # Trás esquerda
]

    wheels = []
    joints = []
    for x, y, z, is_front in wheel_positions:
        wheel_body, wheel_geom = create_wheel(world, space)
        wheel_body.setPosition((x, y, z))
        joint = attach_wheel(world, car_body, wheel_body, x, y, z, is_front)  
        wheels.append((wheel_body, wheel_geom))
        joints.append(joint)
    return car_body, car_geom, wheels, joints

# 🎮 Controlar o carro pelo teclado
def _keyfunc(c, x, y):
    global speed, steer
    if c == b'q': sys.exit(0)
    elif c == b'w': speed += 0.2
    elif c == b's': speed -= 0.2
    elif c == b'a': steer = max(steer - 0.1, -0.5)
    elif c == b'd': steer = min(steer + 0.1, 0.5)
    elif c == b'x': speed, steer = 0, 0

# ⚡ Aplicar movimento ao carro
def update_car():
    for i, joint in enumerate(joints):
        if i < 2:
            joint.setParam(ode.ParamLoStop, steer)
            joint.setParam(ode.ParamHiStop, steer)
        joint.setParam(ode.ParamVel2, speed)
        joint.setParam(ode.ParamFMax2, 1000)

# 🖼️ Renderizar os objetos ODE no OpenGL
def draw_body(body):
    x, y, z = body.getPosition()
    R = body.getRotation()
    rot = [R[0], R[3], R[6], 0., R[1], R[4], R[7], 0., R[2], R[5], R[8], 0., x, y, z, 1.0]
    glPushMatrix()
    glMultMatrixd(rot)
    if body.shape == "box":
        sx, sy, sz = body.boxsize
        glScalef(sx, sy, sz)
        glutSolidCube(1)
    elif body.shape == "cylinder":
        radius, length = body.cyl_params
        glRotatef(90, 0, 1, 0)
        glutSolidCylinder(radius, length, 20, 20)
    glPopMatrix()

# 🎥 Configuração inicial do OpenGL
def prepare_GL():
    glViewport(0, 0, 640, 480)
    glClearColor(0.8, 0.8, 0.9, 0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glEnable(GL_NORMALIZE)
    glShadeModel(GL_FLAT)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45, 1.3333, 0.2, 20000)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(4, 20, 10, 0, 0, 0, 0, 1, 0)

# 🖼️ Callback para desenhar a cena
def _drawfunc():
    prepare_GL()
    draw_body(car_body)
    for wheel, _ in wheels:
        draw_body(wheel)
    glutSwapBuffers()

# 🕹️ Loop de atualização da física e renderização
def _idlefunc():
    global lasttime
    t = dt - (time.time() - lasttime)
    if t > 0: time.sleep(t)
    update_car()
    glutPostRedisplay()
    for _ in range(4):
        space.collide((world, contactgroup), near_callback)
        world.step(dt/4)
        contactgroup.empty()
    lasttime = time.time()

# 💥 Função de colisão
def near_callback(args, geom1, geom2):
    world, contactgroup = args
    contacts = ode.collide(geom1, geom2)
    for c in contacts:
        c.setBounce(0.2)
        c.setMu(5000)
        j = ode.ContactJoint(world, contactgroup, c)
        j.attach(geom1.getBody(), geom2.getBody())

# 🎬 Inicializar GLUT
glutInit([])
glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE)
glutInitWindowSize(640, 480)
glutCreateWindow("Simulação de Carro no ODE")
glutKeyboardFunc(_keyfunc)
glutDisplayFunc(_drawfunc)
glutIdleFunc(_idlefunc)

# 🚗 Criar o carro
car_body, car_geom, wheels, joints = create_car(world, space)

# 🕒 Definições do loop
fps = 50
dt = 1.0 / fps
lasttime = time.time()

# 🎮 Iniciar a simulação
glutMainLoop()

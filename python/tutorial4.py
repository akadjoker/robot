
import sys, os, random, time
from math import *
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

import ode

# geometric utility functions
def scalp (vec, scal):
    vec[0] *= scal
    vec[1] *= scal
    vec[2] *= scal

def length (vec):
    return sqrt (vec[0]**2 + vec[1]**2 + vec[2]**2)

# prepare_GL
def prepare_GL():
    """Prepare drawing.
    """

    # Viewport
    glViewport(0,0,640,480)

    # Initialize
    glClearColor(0.8,0.8,0.9,0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST)
    glDisable(GL_LIGHTING)
    glEnable(GL_LIGHTING)
    glEnable(GL_NORMALIZE)
    glShadeModel(GL_FLAT)

    # Projection
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective (45,1.3333,0.2,20)

    # Initialize ModelView matrix
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

    # Light source
    glLightfv(GL_LIGHT0,GL_POSITION,[0,0,1,0])
    glLightfv(GL_LIGHT0,GL_DIFFUSE,[1,1,1,1])
    glLightfv(GL_LIGHT0,GL_SPECULAR,[1,1,1,1])
    glEnable(GL_LIGHT0)

    # View transformation
    gluLookAt (2.4, 3.6, 4.8, 0.5, 0.5, 0, 0, 1, 0)


def draw_axes():
    glBegin(GL_LINES)
    # Eixo X (vermelho)
    glColor3f(1, 0, 0)
    glVertex3f(0, 0, 0)
    glVertex3f(2, 0, 0)

    # Eixo Y (verde)
    glColor3f(0, 1, 0)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 2, 0)

    # Eixo Z (azul)
    glColor3f(0, 0, 1)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 0, 2)

    glEnd()
    glColor3f(1, 1, 1)  # Restaurar branco

# draw_body
def draw_body(body):
    """Draw an ODE body.
    """

    x,y,z = body.getPosition()
    R = body.getRotation()
    rot = [R[0], R[3], R[6], 0.,
           R[1], R[4], R[7], 0.,
           R[2], R[5], R[8], 0.,
           x, y, z, 1.0]
    glPushMatrix()
    glMultMatrixd(rot)
    if body.shape=="box":
        sx,sy,sz = body.boxsize
        glScalef(sx, sy, sz)
        glutSolidCube(1)
    elif body.shape == "sphere":
        radius = body.radius
        glutSolidSphere(radius, 20, 20)
    elif body.shape == "cylinder":
        radius, length = body.cyl_params
        glutSolidCylinder(radius, length, 20, 20)

    
    draw_axes()
    glPopMatrix()


# create_box
def create_box(world, space, density, lx, ly, lz):
    """Create a box body and its corresponding geom."""

    # Create body
    body = ode.Body(world)
    M = ode.Mass()
    M.setBox(density, lx, ly, lz)
    body.setMass(M)

    # Set parameters for drawing the body
    body.shape = "box"
    body.boxsize = (lx, ly, lz)

    # Create a box geom for collision detection
    geom = ode.GeomBox(space, lengths=body.boxsize)
    geom.setBody(body)

    return body, geom

def create_sphere(world, space, density, radius):
    body = ode.Body(world)
    M = ode.Mass()
    M.setSphere(density, radius)
    body.setMass(M)
    body.shape = "sphere"
    body.radius = radius
    geom = ode.GeomSphere(space, radius)
    geom.setBody(body)
    return body, geom

def create_cylinder(world, space, density, radius, length):
    body = ode.Body(world)
    M = ode.Mass()
    M.setCylinder(density, 3, radius, length)  # 3 = eixo Z
    body.setMass(M)
    body.shape = "cylinder"
    body.cyl_params = (radius, length)
    geom = ode.GeomCylinder(space, radius, length)
    geom.setBody(body)
    return body, geom

# drop_object
def drop_object():
    """Drop an object into the scene."""

    global bodies, geom, counter, objcount

    body, geom = create_box(world, space, 1000, 1.0,0.2,0.2)
    body.setPosition( (random.gauss(0,0.1),3.0,random.gauss(0,0.1)) )
    theta = random.uniform(0,2*pi)
    ct = cos (theta)
    st = sin (theta)
    body.setRotation([ct, 0., -st, 0., 1., 0., st, 0., ct])
    bodies.append(body)
    geoms.append(geom)
    counter=0
    objcount+=1

# explosion
def explosion():
    global bodies

    for b in bodies:
        l=b.getPosition ()
        d = length (l)
        a = max(0, 40000*(1.0-0.2*d*d))
        l = [l[0] / 4, l[1], l[2] /4]
        scalp (l, a / length (l))
        b.addForce(l)

# pull
def pull():
    """Pull the objects back to the origin.

    Every object will be pulled back to the origin.
    Every couple of frames there'll be a thrust upwards so that
    the objects won't stick to the ground all the time.
    """
    global bodies, counter

    for b in bodies:
        l=list (b.getPosition ())
        scalp (l, -1000 / length (l))
        b.addForce(l)
        if counter%60==0:
            b.addForce((0,10000,0))

# Collision callback
def near_callback(args, geom1, geom2):
    # Check if the objects do collide
    contacts = ode.collide(geom1, geom2)

    # Create contact joints
    world,contactgroup = args
    for c in contacts:
        print(f"Colisão entre {geom1} e {geom2}")
      
        c.setBounce(0.2)
        c.setMu(5000)
        j = ode.ContactJoint(world, contactgroup, c)
        j.attach(geom1.getBody(), geom2.getBody())



######################################################################

# Initialize Glut
glutInit ([])

# Open a window
glutInitDisplayMode (GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE)

x = 0
y = 0
width = 640
height = 480
glutInitWindowPosition (x, y);
glutInitWindowSize (width, height);
glutCreateWindow ("testode")

# Create a world object
world = ode.World()
world.setGravity( (0,-9.81,0) )
world.setERP(0.8)
world.setCFM(1E-5)

# Create a space object
space = ode.Space()

# Create a plane geom which prevent the objects from falling forever
floor = ode.GeomPlane(space, (0,1,0), 0)

# A list with ODE bodies
bodies = []

# The geoms for each of the bodies
geoms = []

# A joint group for the contact joints that are generated whenever
# two bodies collide
contactgroup = ode.JointGroup()

# Some variables used inside the simulation loop
fps = 50
dt = 1.0/fps
running = True
state = 0
counter = 0
objcount = 0
lasttime = time.time()


# keyboard callback
# Atualizar a função de teclado
def _keyfunc(c, x, y):
    global objcount

    if c == b'q':  # Sair do programa
        sys.exit(0)

    elif c == b'c':  # Criar caixa
        print("Caixa criada!")
        body, geom = create_box(world, space, 1000, 1.0, 0.2, 0.2)
        body.setPosition((random.gauss(0, 0.1), 3.0, random.gauss(0, 0.1)))
        bodies.append(body)
        geoms.append(geom)
        objcount += 1

    elif c == b's':  # Criar esfera
        print("Esfera criada!")
        body, geom = create_sphere(world, space, 1000, 0.5)
        body.setPosition((random.gauss(0, 0.1), 3.0, random.gauss(0, 0.1)))
        bodies.append(body)
        geoms.append(geom)
        objcount += 1

    elif c == b'y':  # Criar cilindro
        print("Cilindro criado!")
        body, geom = create_cylinder(world, space, 1000, 0.3, 1.0)
        body.setPosition((random.gauss(0, 0.1), 3.0, random.gauss(0, 0.1)))
        bodies.append(body)
        geoms.append(geom)
        objcount += 1

    elif c == b'e':  # Explosão!
        print("Explosão ativada!")
        explosion()

    elif c == b'p':  # Puxar objetos de volta
        print("Puxar objetos para o centro!")
        pull()

glutKeyboardFunc(_keyfunc)


glutKeyboardFunc (_keyfunc)

# draw callback
def _drawfunc ():
    # Draw the scene
    prepare_GL()
    for b in bodies:
        draw_body(b)

    glutSwapBuffers ()

glutDisplayFunc (_drawfunc)

# idle callback
def _idlefunc ():
    global counter, state, lasttime

    t = dt - (time.time() - lasttime)
    if (t > 0):
        time.sleep(t)



    glutPostRedisplay ()

    # Simulate
    n = 4

    for i in range(n):
        # Detect collisions and create contact joints
        space.collide((world,contactgroup), near_callback)

        # Simulation step
        world.step(dt/n)

        # Remove all contact joints
        contactgroup.empty()

    lasttime = time.time()

glutIdleFunc (_idlefunc)

glutMainLoop ()


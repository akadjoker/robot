import ode

# Criar um mundo no ODE
world = ode.World()
world.setGravity((0, -9.81, 0))

# Criar um corpo no mundo
body = ode.Body(world)
mass = ode.Mass()
mass.setSphere(2500, 0.05)  # Massa e raio da esfera
body.setMass(mass)

print("ODE está a funcionar corretamente!")

from dolfin import *

mesh = UnitSquareMesh(3, 3)
V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(v), grad(u)) * dx

A = assemble(a)
from IPython import embed; embed()

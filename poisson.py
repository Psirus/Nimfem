from dolfin import *

mesh = UnitSquareMesh(1000, 1000)
V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(v), grad(u)) * dx

A = assemble(a)
print("Number nonzeros: ", A.nnz())
print("Number elements: ", mesh.num_cells())

from dolfin import *

mesh = UnitSquareMesh(2, 2)
V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(v), grad(u)) * dx
f = Constant(-6.0)
L = f*v*dx

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

A = assemble(a)
rhs = assemble(L)
print(A.array())
print(rhs.get_local())
# Compute solution
u = Function(V)
solve(a == L, u, bc)

# A = assemble(a)
vtkfile = File('fenics_output.pvd')
vtkfile << u

print(mesh.num_cells())

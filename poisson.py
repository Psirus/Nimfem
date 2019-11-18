import sys
from dolfin import *

mesh = UnitSquareMesh(50, 50)
V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(v), grad(u)) * dx
f = Constant(-6.0)
L = f * v * dx

A = assemble(a)
sys.exit(0)

# Define boundary condition
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)


def boundary(x, on_boundary):
    return on_boundary


bc = DirichletBC(V, u_D, boundary)

u = Function(V)
solve(
    a == L,
    u,
    bc,
    solver_parameters={"linear_solver": "cg"},
    form_compiler_parameters={"quadrature_degree": 2},
)

vtkfile = File("fenics_output.pvd")
vtkfile << u

print("Number of elements: ", mesh.num_cells())

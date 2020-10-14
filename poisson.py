import sys
from dolfin import *

mesh = UnitSquareMesh(50, 50)
V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(v), grad(u)) * dx
f = Constant(-6.0)
L = f * v * dx

# Define boundary condition
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)


def boundary(x, on_boundary):
    return on_boundary


bc = DirichletBC(V, u_D, boundary)

u = Function(V)
problem = LinearVariationalProblem(a, L, u, bc, form_compiler_parameters={"quadrature_degree": 2})
solver = LinearVariationalSolver(problem)
prms = solver.parameters
prms["krylov_solver"]["monitor_convergence"] = False
prms["krylov_solver"]["relative_tolerance"] = 1e-16
prms["krylov_solver"]["absolute_tolerance"] = 1e-9
prms["linear_solver"] = "cg"
prms["preconditioner"] = "ilu"

solver.solve()

vtkfile = File("fenics_output.pvd")
vtkfile << u

print("Number of elements: ", mesh.num_cells())

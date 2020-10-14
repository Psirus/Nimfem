import mesh
import triangle
import dense
import assembly
import sparse
import iterative_methods
import io
import times

proc myF(x: Vector[2], J: Matrix[2, 2]): Matrix[3, 3] =
  let B = inv(J) * derivativeShapeFunctions()
  result = B.transpose * B

proc mySource(x: Vector[2], J: Matrix[2, 2]): Vector[3] =
  let N = shapeFunctions(x)
  result = - 6.0 * N

proc bc(x: Vector[2]): float =
  result = 1.0 + x[0]*x[0] + 2.0*x[1]*x[1]

let my_mesh = UnitSquareMesh(50)

var A = assembleMatrix(myF, my_mesh)
var f = assembleVector(mySource, my_mesh)

setDiagonalRows(A, my_mesh.boundary_nodes)
applyBC(f, my_mesh, bc)

var before = cpuTime()
var u = conjugate_gradient(A, f)
echo "CG time:", cpuTime() - before

before = cpuTime()
var P = incomplete_lu(A)
echo "ILU time:", cpuTime() - before

before = cpuTime()
var u2 = preconditioned_cg(A, P, f)
echo "PCG time:", cpuTime() - before

writeVTK(my_mesh, u)

echo "Number of elements: ", my_mesh.connectivity.len

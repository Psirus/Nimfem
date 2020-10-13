import mesh
import triangle
import dense
import assembly
import sparse
import iterative_methods
import io

proc myF(x: Vector[2], J: Matrix[2, 2]): Matrix[3, 3] =
  let B = inv(J) * derivativeShapeFunctions()
  result = B.transpose * B

proc mySource(x: Vector[2], J: Matrix[2, 2]): Vector[3] =
  let N = shapeFunctions(x)
  result = - 6.0 * N

proc bc(x: Vector[2]): float =
  result = 1.0 + x[0]*x[0] + 2.0*x[1]*x[1]

let my_mesh = UnitSquareMesh(3)

var A = assembleMatrix(myF, my_mesh)
var f = assembleVector(mySource, my_mesh)

setDiagonalRows(A, my_mesh.boundary_nodes)
applyBC(f, my_mesh, bc)

var u = conjugate_gradient(A, f)
var P = incomplete_lu2(A)

P.setEntry(    6,   5,  -0.25000)
P.setEntry(   10,   6,  -0.26667)
P.setEntry(   9,  9,   3.75000)
P.setEntry(   10,  9,  -0.26667)
P.setEntry(   10,  10,   3.46667)

var u2 = preconditioned_cg(A, P, f)

writeVTK(my_mesh, u)

echo "Number of elements: ", my_mesh.connectivity.len

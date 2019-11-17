import mesh
import triangle
import dense
import assembly
import sparse
import iterative_methods
import io

proc myF(x: Vector, J: Matrix): Matrix =
  let B = inv(J) * derivativeShapeFunctions()
  result = B.transpose * B

proc mySource(x: Vector, J: Matrix): Vector =
  let N = shapeFunctions(x)
  result = - 6.0 * N

proc bc(x: Vector): float =
  result = 1.0 + x[0]*x[0] + 2.0*x[1]*x[1]

let my_mesh = UnitSquareMesh(700)
echo my_mesh.connectivity.len

var A = assembleMatrix(myF, my_mesh)
var f = assembleVector(mySource, my_mesh)

setDiagonalRows(A, my_mesh.boundary_nodes)
applyBC(f, my_mesh, bc)

var u = conjugate_gradient(A, f)

writeVTK(my_mesh, u)

echo "Number nonzeros: ", A.aa.len
echo "Number elements: ", my_mesh.connectivity.len

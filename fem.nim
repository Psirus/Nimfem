import mesh
import triangle
import dense
import assembly
import sparse

proc myF(x: Vector): Matrix =
  let B = derivativeShapeFunctions()
  result = B.transpose * B

proc mySource(x: Vector): Vector =
  let N = shapeFunctions(x)
  result = 6 * N

proc bc(x: Vector): float =
  result = 1.0 + x[0]*x[0] + 2.0*x[1]*x[1]

let my_mesh = UnitSquareMesh(3)
var A = assembleMatrix(myF, my_mesh)
echo "Number nonzeros: ", A.aa.len
echo "Number elements: ", my_mesh.connectivity.len

var f = assembleVector(mySource, my_mesh)
echo f

setDiagonalRows(A, my_mesh.boundary_nodes)
applyBC(f, my_mesh, bc)

echo f

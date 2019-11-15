import dense
import sparse
import mesh
import quadrature
import triangle

proc myF(x: Vector): Matrix =
  let B = derivativeShapeFunctions()
  result = B.transpose * B

proc assemble(f: proc(x: Vector): Matrix, mesh: Mesh): SparseMatrix =
  var triplets = newSeq[(int, int, float)]()
  for elem in mesh.connectivity:
    var localNodes = newSeqOfCap[array[2, float]](3)
    for dof in elem:
      localNodes.add(mesh.nodes[dof])
    let elem_matrix = triangleSecondOrderQuadrature(f, localNodes)
    for i, dof_i in elem.pairs:
      for j, dof_j in elem.pairs:
        triplets.add((dof_i, dof_j, elem_matrix[i, j]))

  result = fromTripletList(triplets)
    

let my_mesh = UnitSquareMesh(3)
let A = assemble(myF, my_mesh)
echo A.ia
echo A.ja
echo A.aa

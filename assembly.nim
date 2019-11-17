import sequtils

import dense
import sparse
import mesh
import quadrature

proc assembleMatrix*(f: proc(x: Vector, J: Matrix): Matrix, mesh: Mesh): SparseMatrix =
  let num_nodes = mesh.nodes.len
  var Ai = newSeqOfCap[int](15*num_nodes)
  var Aj = newSeqOfCap[int](15*num_nodes)
  var Ax = newSeqOfCap[float](15*num_nodes)
  for elem in mesh.connectivity:
    var localNodes = newSeqOfCap[array[2, float]](3)
    for dof in elem:
      localNodes.add(mesh.nodes[dof])
    let elem_matrix = triangleSecondOrderQuadrature(f, localNodes)
    for i, dof_i in elem.pairs:
      for j, dof_j in elem.pairs:
        Ai.add(dof_i)
        Aj.add(dof_j)
        Ax.add(elem_matrix[i, j])

  result = toCSR(Ai, Aj, Ax)

proc assembleVector*(f: proc(x: Vector, J: Matrix): Vector, mesh: Mesh): Vector =
  result = newVector(mesh.nodes.len)
  for elem in mesh.connectivity:
    var localNodes = newSeqOfCap[array[2, float]](3)
    for dof in elem:
      localNodes.add(mesh.nodes[dof])
    let elem_vec = triangleSecondOrderQuadratureVector(f, localNodes)
    for i, dof_i in elem.pairs:
      result[dof_i] += elem_vec[i]

proc applyBC*(f: var Vector, mesh: Mesh, bc: proc(x: Vector): float) =
  for node in mesh.boundary_nodes:
    let coordinates = toSeq(mesh.nodes[node])
    f[node] = bc(coordinates)

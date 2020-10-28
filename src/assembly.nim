import dense
import sparse
import mesh
import quadrature

proc assembleMatrix*[M, N, D](f: proc(x: Vector[D], J: Matrix[D, D]):
    Matrix[M, N], mesh: Mesh): SparseMatrix =
  ## Assemble element matrices into global matrix.
  ## The element matrices are the result of integrating `f` over the elemnt.
  let num_nodes = mesh.nodes.len
  var Ai = newSeqOfCap[int](15*num_nodes)
  var Aj = newSeqOfCap[int](15*num_nodes)
  var Ax = newSeqOfCap[float](15*num_nodes)
  for elem in mesh.connectivity:
    var localNodes: array[3, array[2, float]]
    for i, dof in elem.pairs:
      localNodes[i] = mesh.nodes[dof]
    let elem_matrix = triSecondOrderQuadrature(f, localNodes)
    for i, dof_i in elem.pairs:
      for j, dof_j in elem.pairs:
        Ai.add(dof_i)
        Aj.add(dof_j)
        Ax.add(elem_matrix[i, j])

  result = toCSR(Ai, Aj, Ax)

proc assembleVector*[N, D](f: proc(x: Vector[D], J: Matrix[D, D]): Vector[N],
    mesh: Mesh): DynamicVector =
  ## Assemble element vectors into global vector.
  ## The element matrices are the result of integrating `f` over the elemnt.
  result = newVector(mesh.nodes.len)
  for elem in mesh.connectivity:
    var localNodes: array[N, array[D, float]]
    for i, dof in elem.pairs:
      localNodes[i] = mesh.nodes[dof]
    let elem_vec = triSecondOrderQuadratureVec(f, localNodes)
    for i, dof_i in elem.pairs:
      result[dof_i] += elem_vec[i]

proc applyBC*(f: var DynamicVector, mesh: Mesh, bc: proc(x: Vector[2]): float) =
  ## Apply the boundary condition `bc` to the global vector `f` at all boundary nodes.
  for node in mesh.boundary_nodes:
    let coordinates = mesh.nodes[node]
    f[node] = bc(coordinates)

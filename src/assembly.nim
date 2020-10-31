import dense
import sparse
import mesh
import quadrature

proc assembleMatrix*[N, D, Element](f: proc(x: Vector[D], J: Matrix[D, D]):
    Matrix[N, N], mesh: Mesh[Element]): SparseMatrix =
  ## Assemble element matrices into global matrix.
  ## The element matrices are the result of integrating `f` over the elemnt.
  let num_dofs=3
  let num_entries = N*N*mesh.connectivity.len*num_dofs
  var Ai = newSeqOfCap[int](num_entries)
  var Aj = newSeqOfCap[int](num_entries)
  var Ax = newSeqOfCap[float](num_entries)
  for elem in mesh.connectivity:
    var localNodes: array[Element.num_nodes(), array[D, float]]
    for i, dof in elem.pairs:
      localNodes[i] = mesh.nodes[dof]
    when Element.shape == triangle:
      let elem_matrix = triSecondOrderQuadrature(f, localNodes)
    elif Element.shape == tetrahedron:
      let elem_matrix = tetSecondOrderQuadrature(f, localNodes)
    #echo dense.`$`(elem_matrix)
    for i, dof_i in elem.pairs:
      for j, dof_j in elem.pairs:
        for ii in 0..<num_dofs:
          for jj in 0..<num_dofs:
            Ai.add(num_dofs*dof_i+ii)
            Aj.add(num_dofs*dof_j+jj)
            Ax.add(elem_matrix[i+ii, j+jj])

  result = toCSR(Ai, Aj, Ax)

proc assembleVector*[N, D, Element](f: proc(x: Vector[D], J: Matrix[D, D]): Vector[N],
    mesh: Mesh[Element]): DynamicVector =
  ## Assemble element vectors into global vector.
  ## The element matrices are the result of integrating `f` over the elemnt.
  let num_dofs=3
  result = newVector(num_dofs*mesh.nodes.len)
  for elem in mesh.connectivity:
    var localNodes: array[Element.num_nodes(), array[D, float]]
    for i, dof in elem.pairs:
      localNodes[i] = mesh.nodes[dof]
    when Element.shape == triangle:
      let elem_vec = triSecondOrderQuadratureVec(f, localNodes)
    elif Element.shape == tetrahedron:
      let elem_vec = tetSecondOrderQuadratureVec(f, localNodes)
    for i, dof_i in elem.pairs:
      for ii in 0..<num_dofs:
        result[num_dofs*dof_i + ii] += elem_vec[num_dofs*i+ii]

proc applyBC*(f: var DynamicVector, mesh: Mesh, bc: proc(x: Vector[3]): Vector[3], boundary_nodes: seq[int]) =
  ## Apply the boundary condition `bc` to the global vector `f` at all boundary nodes.
  let num_dofs = 3
  for node in boundary_nodes:
    let coordinates = mesh.nodes[node]
    let bcs = bc(coordinates)
    for ii in 0..<num_dofs:
      f[num_dofs*node + ii] = bcs[ii]

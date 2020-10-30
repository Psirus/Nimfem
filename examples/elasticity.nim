import ../nimfem

proc elasticityTensor(E, nu: float): Matrix[6, 6] =
  let factor = E / ((1.0 + nu) * (1.0 - 2.0*nu))
  result = [[1.0 - nu, nu, nu, 0.0, 0.0, 0.0],
            [nu, 1.0 - nu, nu, 0.0, 0.0, 0.0],
            [nu, nu, 1.0 - nu, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0 - nu, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0 - nu, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0 - nu]]
  result = result * factor

proc myF(x: Vector[3], J: Matrix[3, 3]): Matrix[12, 12] =
  let
   mu = 1.0
   lambda = 1.25
   E = mu * (3.0*lambda + 2.0*mu) / (lambda + mu)
   nu = lambda + 2.0 * mu

   C = elasticityTensor(E, nu)
   invJ = inv(J)
   B = kelvinStrain(invJ)

  result = B.transpose * C * B

proc mySource(x: Vector[3], J: Matrix[3, 3]): Vector[12] =
  let N = vectorShapeFunctions(x)
  let rho = 1.0
  let g = 2.0 / 125.0
  let f = [0.0, 0.0, -rho * g]
  result = transpose(N) * f

proc findLeftBoundary(mesh: Mesh): seq[int] =
  for i, node in mesh.nodes.pairs:
    if node[0] == 0.0:
      result.add(i)

proc bc(x: Vector[3]): Vector[3] =
  result = [0.0, 0.0, 0.0]

let my_mesh = readDolfinXml(LagrangeTetrahedron(1), "test/data/beam.xml")

let boundary_nodes = findLeftBoundary(my_mesh)

var A = assembleMatrix(myF, my_mesh)
var f = assembleVector(mySource, my_mesh)

setDiagonalRows(A, boundary_nodes)
applyBC(f, my_mesh, bc, boundary_nodes)

var P = incomplete_lu(A)
var u = preconditioned_cg(A, P, f)

writeVTK("output.vtk", my_mesh, u)

echo "Number of elements: ", my_mesh.connectivity.len

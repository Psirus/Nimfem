import math, xmlparser, xmltree, strutils

type
  Element = object
    num_nodes: int
    dim: int

proc num_nodes*(e: Element): int = e.num_nodes
proc dim*(e: Element): int = e.dim

type
  Mesh*[element: static[Element]] = object
    nodes*: seq[array[element.dim, float]]
    connectivity*: seq[array[element.num_nodes, int]]
    boundary_nodes*: seq[int]

proc LagrangeLine*(order: int): Element =
  result.num_nodes = order + 1
  result.dim = 1

proc LagrangeTriangle*(order: int): Element =
  result.num_nodes = (order + 2) * (order + 1) div 2
  result.dim = 2

proc LagrangeTetrahedron*(order: int): Element =
  result.num_nodes = (order + 3) * (order + 2) * (order + 1) div 6
  result.dim = 3

proc UnitSquareMesh*(n: int): Mesh[LagrangeTriangle(1)] =
  ## Create a triangle mesh of a unit square.
  ## n is the number of elements per edge.
  result.nodes = newSeq[array[2, float]]((n + 1)^2)
  for i in 0..n:
    for j in 0..n:
      if (i == 0) or (j == 0) or (i == n) or (j == n):
        result.boundary_nodes.add(i + j*(n+1))
      let x = float(i) / float(n)
      let y = 1.0 - float(j) / float(n)
      result.nodes[j * (n + 1) + i] = [x, y]

  result.connectivity = newSeqOfCap[array[3, int]](2 * n^2)
  for j in 0..<n:
    for i in 0..<n:
      result.connectivity.add([i + j*(n+1), i + (j+1)*(n+1), i + (j+1)*(n+1) + 1])
      result.connectivity.add([i + j*(n+1), i + (j+1)*(n+1) + 1, i + j*(n+1) + 1])

proc UnitIntervalMesh*(n: int): Mesh[LagrangeLine(1)] =
  ## Create a mesh of a unit interval.
  ## n is the number of elements.
  result.nodes = newSeq[array[1, float]](n+1)
  for i in 0..n:
    let x = float(i) / float(n)
    result.nodes[i] = [x]

  result.boundary_nodes.add(0)
  result.boundary_nodes.add(n)

  result.connectivity = newSeqOfCap[array[2, int]](n)
  for i in 0..<n:
    result.connectivity.add([i, i+1])

proc transform*(mesh: Mesh, f: proc (x: float): float): Mesh =
  ## Transform a 1D mesh using a function `f`.
  result = mesh
  for node in result.nodes.mItems():
    for i in 0..<node.len:
      node[i] = f(node[i])

proc readDolfinXml*(element: static[Element], filename: string): Mesh[element] =
  ## Read in Dolfin XML meshes from `filename`.
  ## For compile-time optimizations, the elementtype of the resulting mesh is required as `element`.
  let xml_contents = loadXml(filename)
  let mesh_element = child(xml_contents, "mesh")

  let accetable_celltypes = ["triangle", "tetrahedron"]
  if mesh_element.attr("celltype") notin accetable_celltypes:
    quit("Celltype currently not supported.")

  let vertices_element = child(mesh_element, "vertices")
  let num_nodes = parseInt(vertices_element.attr("size"))
  result.nodes = newSeq[array[element.dim, float]](num_nodes)
  for vertex in vertices_element.items():
    let index = parseInt(vertex.attr("index"))
    let x = parseFloat(vertex.attr("x"))
    let y = parseFloat(vertex.attr("y"))
    when element.dim == 2:
      result.nodes[index] = [x, y]
    elif element.dim == 3:
      let z = parseFloat(vertex.attr("z"))
      result.nodes[index] = [x, y, z]

  let cells_element = child(mesh_element, "cells")
  let num_elements = parseInt(cells_element.attr("size"))
  result.connectivity = newSeq[array[element.num_nodes, int]](num_elements)
  for cell in cells_element.items():
    let index = parseInt(cell.attr("index"))
    let v0 = parseInt(cell.attr("v0"))
    let v1 = parseInt(cell.attr("v1"))
    let v2 = parseInt(cell.attr("v2"))
    when element.num_nodes == 3:
      result.connectivity[index] = [v0, v1, v2]
    elif element.num_nodes == 4:
      let v3 = parseInt(cell.attr("v3"))
      result.connectivity[index] = [v0, v1, v2, v3]

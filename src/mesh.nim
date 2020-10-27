import math, xmlparser, xmltree, strutils

type
  Mesh*[Dim, NodesPerElement: static[int]]  = object
    nodes*: seq[array[Dim, float]]
    connectivity*: seq[array[NodesPerElement, int]]
    boundary_nodes*: seq[int]

# division means number of elements per edge of the unit square
proc UnitSquareMesh*(n: int): Mesh[2, 3] =
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


proc readDolfinXml*(D, N: static[int], filename: string): Mesh[D, N] =
  let xml_contents = loadXml(filename)
  let mesh_element = child(xml_contents, "mesh")

  let accetable_celltypes = ["triangle", "tetrahedron"]
  if mesh_element.attr("celltype") notin accetable_celltypes:
    quit("Celltype currently not supported.")

  let vertices_element = child(mesh_element, "vertices")
  let num_nodes = parseInt(vertices_element.attr("size"))
  result.nodes = newSeq[array[D, float]](num_nodes)
  for vertex in vertices_element.items():
    let index = parseInt(vertex.attr("index"))
    let x = parseFloat(vertex.attr("x"))
    let y = parseFloat(vertex.attr("y"))
    when D == 2:
      result.nodes[index] = [x, y]
    elif D == 3:
      let z = parseFloat(vertex.attr("z"))
      result.nodes[index] = [x, y, z]

  let cells_element = child(mesh_element, "cells")
  let num_elements = parseInt(cells_element.attr("size"))
  result.connectivity = newSeq[array[N, int]](num_elements)
  for cell in cells_element.items():
    let index = parseInt(cell.attr("index"))
    let v0 = parseInt(cell.attr("v0"))
    let v1 = parseInt(cell.attr("v1"))
    let v2 = parseInt(cell.attr("v2"))
    when N == 3:
      result.connectivity[index] = [v0, v1, v2]
    elif N == 4:
      let v3 = parseInt(cell.attr("v3"))
      result.connectivity[index] = [v0, v1, v2, v3]

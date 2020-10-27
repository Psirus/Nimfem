import math, xmlparser, xmltree, strutils

type
  Mesh* = object
    nodes*: seq[array[2, float]]
    connectivity*: seq[array[3, int]]
    boundary_nodes*: seq[int]

# division means number of elements per edge of the unit square
proc UnitSquareMesh*(n: int): Mesh =
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

proc readDolfinXml*(filename: string): Mesh =
  let xml_contents = loadXml(filename)
  let mesh_element = child(xml_contents, "mesh")
  if mesh_element.attr("celltype") != "triangle":
    quit("Celltype currently not supported.")

  let vertices_element = child(mesh_element, "vertices")
  let num_nodes = parseInt(vertices_element.attr("size"))
  result.nodes = newSeq[array[2, float]](num_nodes)
  for vertex in vertices_element.items():
    let index = parseInt(vertex.attr("index"))
    let x = parseFloat(vertex.attr("x"))
    let y = parseFloat(vertex.attr("y"))
    result.nodes[index] = [x, y]

  let cells_element = child(mesh_element, "cells")
  let num_elements = parseInt(cells_element.attr("size"))
  result.connectivity = newSeq[array[3, int]](num_elements)
  for cell in cells_element.items():
    let index = parseInt(cell.attr("index"))
    let v0 = parseInt(cell.attr("v0"))
    let v1 = parseInt(cell.attr("v1"))
    let v2 = parseInt(cell.attr("v2"))
    result.connectivity[index] = [v0, v1, v2]

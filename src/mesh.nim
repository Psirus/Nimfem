import math

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


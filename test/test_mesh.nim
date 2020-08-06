import ../mesh

block:
  let mesh = UnitSquareMesh(3)
  doAssert mesh.nodes.len == 16
  doAssert mesh.connectivity.len == 18

block:
  let mesh = UnitSquareMesh(1)
  doAssert mesh.nodes == @[[0.0, 1.0], [1.0, 1.0], [0.0, 0.0], [1.0, 0.0]]
  doAssert mesh.connectivity == @[[0, 2, 3], [0, 3, 1]]

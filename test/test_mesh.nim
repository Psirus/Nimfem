import ../src/mesh

block:
  let mesh = UnitSquareMesh(3)
  doAssert mesh.nodes.len == 16
  doAssert mesh.connectivity.len == 18

block:
  let mesh = UnitSquareMesh(1)
  doAssert mesh.nodes == @[[0.0, 1.0], [1.0, 1.0], [0.0, 0.0], [1.0, 0.0]]
  doAssert mesh.connectivity == @[[0, 2, 3], [0, 3, 1]]

block:
  let mesh = readDolfinXml("data/ell_2d.xml")

  doAssert mesh.nodes.len == 26
  doAssert mesh.nodes[17][0] == 0.49609
  doAssert mesh.nodes[17][1] == 1.29651

  doAssert mesh.connectivity.len == 34
  doAssert mesh.connectivity[7] == [8, 18, 7]

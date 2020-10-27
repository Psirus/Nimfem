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
  let mesh = readDolfinXml(2, 3, "test/data/ell_2d.xml")

  doAssert mesh.nodes.len == 26
  doAssert mesh.nodes[17][0] == 0.49609
  doAssert mesh.nodes[17][1] == 1.29651

  doAssert mesh.connectivity.len == 34
  doAssert mesh.connectivity[7] == [8, 18, 7]

block:
  let mesh = readDolfinXml(3, 4, "test/data/beam.xml")

  doAssert mesh.nodes.len == 176
  doAssert mesh.nodes[17][0] == 6.0e-1
  doAssert mesh.nodes[17][1] == 2.0e-1/3.0

  doAssert mesh.connectivity.len == 540
  doAssert mesh.connectivity[7] == [1, 2, 46, 57]

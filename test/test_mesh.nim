import unittest

import ../src/mesh

suite "Mesh":

  test "Small unit square":
    let mesh = UnitSquareMesh(3)
    check mesh.nodes.len == 16
    check mesh.connectivity.len == 18

  test "Two triangles square":
    let mesh = UnitSquareMesh(1)
    check mesh.nodes == @[[0.0, 1.0], [1.0, 1.0], [0.0, 0.0], [1.0, 0.0]]
    check mesh.connectivity == @[[0, 2, 3], [0, 3, 1]]

  test "Small interval mesh":
    let mesh = UnitIntervalMesh(5)
    check mesh.nodes.len == 6
    check mesh.connectivity.len == 5

    let mesh2 = transform(mesh, proc (x: float): float = 2*x)
    check mesh2.nodes == @[[0.0], [0.4], [0.8], [1.2], [1.6], [2.0]]

  test "Triangle mesh from file":
    let mesh = readDolfinXml(LagrangeTriangle(1), "test/data/ell_2d.xml")

    check mesh.nodes.len == 26
    check mesh.nodes[17][0] == 0.49609
    check mesh.nodes[17][1] == 1.29651

    check mesh.connectivity.len == 34
    check mesh.connectivity[7] == [8, 18, 7]

  test "Tetrahedron mesh from file":
    let mesh = readDolfinXml(LagrangeTetrahedron(1), "test/data/beam.xml")

    check mesh.nodes.len == 176
    check mesh.nodes[17][0] == 6.0e-1
    check mesh.nodes[17][1] == 2.0e-1/3.0

    check mesh.connectivity.len == 540
    check mesh.connectivity[7] == [1, 2, 46, 57]

  test "Linear triangle":
    let element = LagrangeTriangle(1)
    check element.num_nodes == 3
    check element.dim == 2

  test "Quadratic triangle":
    let element = LagrangeTriangle(2)
    check element.num_nodes == 6
    check element.dim == 2

  test "Linear tetrahedron":
    let element = LagrangeTetrahedron(1)
    check element.num_nodes == 4
    check element.dim == 3

  test "Quadratic tetrahedron":
    let element = LagrangeTetrahedron(2)
    check element.num_nodes == 10
    check element.dim == 3

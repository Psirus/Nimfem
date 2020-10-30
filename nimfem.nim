## .. include:: README.rst

import src/dense
export Vector, Matrix, inv, `*`, transpose

import src/mesh
export Mesh, UnitSquareMesh, UnitIntervalMesh, readDolfinXml, LagrangeLine, LagrangeTriangle, LagrangeTetrahedron

import src/triangle
export shapeFunctions, derivativeShapeFunctions

import src/tetrahedron
export shapeFunctions, kelvinStrain, vectorShapeFunctions

import src/assembly
export assembleVector, assembleMatrix, applyBC

import src/sparse
export setDiagonalRows, rows, cols, `$`, getEntry

import src/iterative_methods
export incomplete_lu, preconditioned_cg, conjugate_gradient

import src/io
export writeVTK

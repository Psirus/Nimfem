## .. include:: README.rst

import src/dense
export Vector, Matrix, inv, `*`, transpose

import src/mesh
export UnitSquareMesh

import src/triangle
export shapeFunctions, derivativeShapeFunctions

import src/assembly
export assembleVector, assembleMatrix, applyBC

import src/sparse
export setDiagonalRows

import src/iterative_methods
export incomplete_lu, preconditioned_cg, conjugate_gradient

import src/io
export writeVTK

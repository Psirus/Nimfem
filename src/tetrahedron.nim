import dense
import math

proc shapeFunctions*(x: Vector[3]): Vector[4] =
  ## Linear shapefunctions.
  assert x[0] <= 1.0
  assert x[1] <= 1.0 - x[0]
  assert x[2] <= 1.0 - x[0] - x[1]
  result = [1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2]]

proc derivativeShapeFunctions*(): Matrix[3, 4] =
  ## First derivative of the shapefunctions.
  result = [[-1.0, 1.0, 0.0, 0.0],
            [-1.0, 0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0, 1.0]]

proc jacobian*(nodes: array[4, array[3, float]]): Matrix[3, 3] =
  ## Jacobian matrix of tetrahedron element.
  result = [[nodes[1,0] - nodes[0,0], nodes[1,1] - nodes[0,1], nodes[1,2] - nodes[0,2]],
            [nodes[2,0] - nodes[0,0], nodes[2,1] - nodes[0,1], nodes[2,2] - nodes[0,2]],
            [nodes[3,0] - nodes[0,0], nodes[3,1] - nodes[0,1], nodes[3,2] - nodes[0,2]]]


proc vectorShapeFunctions*(x: Vector[3]): Matrix[3, 12] =
  let scalarShapeFunctions = shapeFunctions(x)
  for i in 0 .. 2:
    result[0, 3*i] = scalarShapeFunctions[i]
    result[1, 3*i+1] = scalarShapeFunctions[i]
    result[2, 3*i+2] = scalarShapeFunctions[i]
     
# should maybe be named gradient or some such
proc kelvinStrain*(invJ: Matrix[3, 3]): Matrix[6, 12] =
  let a = 1.0 / sqrt(2.0)
  let scalarB = invJ * derivativeShapeFunctions()
  for i in 0..3:
    # σ_11, σ_22, σ_33 (or ε_ii)
    result[0, 3*i]   = scalarB[0, i]
    result[1, 3*i+1] = scalarB[1, i]
    result[2, 3*i+2] = scalarB[2, i]

    # σ_23, σ_13, σ_12 (or ε_ij)
    result[3, 3*i+1] = a * scalarB[2, i]
    result[3, 3*i+2] = a * scalarB[1, i]

    result[4, 3*i]   = a * scalarB[2, i]
    result[4, 3*i+2] = a * scalarB[0, i]

    result[5, 3*i]   = a * scalarB[1, i]
    result[5, 3*i+1] = a * scalarB[0, i]

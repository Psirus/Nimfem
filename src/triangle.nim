import dense

proc shapeFunctions*(x: Vector[2]): Vector[3] =
  ## Linear shapefunctions.
  assert x[0] <= 1.0
  assert x[1] <= 1.0 - x[0]
  result = [1.0 - x[0] - x[1], x[0], x[1]]

proc derivativeShapeFunctions*(): Matrix[2, 3] =
  ## First derivative of the shapefunctions.
  result = [[-1.0, 1.0, 0.0],
            [-1.0, 0.0, 1.0]]

proc jacobian*(nodes: array[3, array[2, float]]): Matrix[2, 2] =
  ## Jacobian matrix of triangle element.
  result = [[nodes[1][0] - nodes[0][0], nodes[1][1] - nodes[0][1]],
            [nodes[2][0] - nodes[0][0], nodes[2][1] - nodes[0][1]]]

proc vectorShapeFunctions*[N: static[int]](x: Vector[2]): Matrix[N, 3*N] =
  let scalarShapeFunctions = shapeFunctions(x)
  for i in 0 ..< N:
    result[0, 3*i] = scalarShapeFunctions[i]
    result[1, 3*i+1] = scalarShapeFunctions[i]
    result[2, 3*i+2] = scalarShapeFunctions[i]
     
# should maybe be named gradient or some such
proc derivativeVectorShapeFunctions*[N: static[int]](invJ: Matrix[2, 2]): Matrix[2*N, 3*N] =
  let scalarB = invJ * derivativeShapeFunctions()
  for i in 0 ..< N:
    result[i, i] = scalarB[0, 0]
    result[i, i+N] = scalarB[0, 1]
    result[i, i+2*N] = scalarB[0, 2]

    result[i + N, i] = scalarB[1, 0]
    result[i + N, i+N] = scalarB[1, 1]
    result[i + N, i+2*N] = scalarB[1, 2]

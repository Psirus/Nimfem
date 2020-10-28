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

import dense

proc shapeFunctions*(x: Vector): Vector =
  assert x[0] <= 1.0
  assert x[1] <= 1.0 - x[0]
  result = @[1.0 - x[0] - x[1], x[0], x[1]]

proc derivativeShapeFunctions*(): Matrix =
  result = @[@[-1.0, 1.0, 0.0],
             @[-1.0, 0.0, 1.0]]

proc jacobian*(nodes: seq[array[2, float]]): Matrix =
  result = @[@[nodes[1][0] - nodes[0][0], nodes[1][1] - nodes[0][1]],
             @[nodes[2][0] - nodes[0][0], nodes[2][1] - nodes[0][1]]]

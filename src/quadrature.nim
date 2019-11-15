import sequtils
import dense
import triangle

# TODO: try sugar syntax
proc triangleSecondOrderQuadrature*(f: proc(x: Vector): Matrix, nodes: seq[array[2, float]]): Matrix =
  let weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
  let gauss_points = [@[1.0/6.0, 1.0/6.0], @[2.0/3.0, 1.0/6.0], @[1.0/6.0, 2.0/3.0]]
  let result_shape = f(@[0.0, 0.0]).shape
  let detJ = determinant(jacobian(nodes))
  result = newMatrix(result_shape)
  for (w, x) in zip(weights, gauss_points):
    result = result + f(x) * (w * 0.5) * detJ

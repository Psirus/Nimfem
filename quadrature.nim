import sequtils
import dense
import triangle

# TODO: try sugar syntax
proc triangleSecondOrderQuadrature*(f: proc(x: Vector, J: Matrix): Matrix, nodes: seq[array[2, float]]): Matrix =
  let weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
  let gauss_points = [@[1.0/6.0, 1.0/6.0], @[2.0/3.0, 1.0/6.0], @[1.0/6.0, 2.0/3.0]]
  let J = jacobian(nodes)
  let detJ = determinant(J)
  let result_shape = f(@[0.0, 0.0], J).shape
  result = newMatrix(result_shape)
  for (w, x) in zip(weights, gauss_points):
    result = result + f(x, J) * (w * 0.5) * detJ

proc triangleSecondOrderQuadratureVector*(f: proc(x: Vector, J: Matrix): Vector, nodes: seq[array[2, float]]): Vector =
  let weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
  let gauss_points = [@[1.0/6.0, 1.0/6.0], @[2.0/3.0, 1.0/6.0], @[1.0/6.0, 2.0/3.0]]
  let J = jacobian(nodes)
  let detJ = determinant(J)
  let result_shape = f(@[0.0, 0.0], J).size
  result = newVector(result_shape)
  for (w, x) in zip(weights, gauss_points):
    result = result + (w * 0.5 * detJ) * f(x, J)

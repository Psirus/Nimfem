import sequtils
import dense
import triangle

# TODO: try sugar syntax
proc triSecondOrderQuadrature*[N, D](f: proc(x: Vector[D], J: Matrix[D, D]):
    Matrix[N, N], nodes: array[N, array[D, float]]): Matrix[N, N] =
  let weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
  let gauss_points = [[1.0/6.0, 1.0/6.0], [2.0/3.0, 1.0/6.0], [1.0/6.0, 2.0/3.0]]
  let J = jacobian(nodes)
  let detJ = determinant(J)
  for (w, x) in zip(weights, gauss_points):
    result = result + f(x, J) * (w * 0.5) * detJ

proc triSecondOrderQuadratureVec*[N, D](f: proc(x: Vector[D], J: Matrix[D, D]):
    Vector[N], nodes: array[N, array[D, float]]): Vector[N] =
  let weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
  let gauss_points = [[1.0/6.0, 1.0/6.0], [2.0/3.0, 1.0/6.0], [1.0/6.0, 2.0/3.0]]
  let J = jacobian(nodes)
  let detJ = determinant(J)
  for (w, x) in zip(weights, gauss_points):
    result = result + (w * 0.5 * detJ) * f(x, J)

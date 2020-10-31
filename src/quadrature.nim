import sequtils
import dense
import triangle
import tetrahedron

# TODO: try sugar syntax
proc triSecondOrderQuadrature*[N, D](f: proc(x: Vector[D], J: Matrix[D, D]):
    Matrix[N, N], nodes: array[3, array[D, float]]): Matrix[N, N] =
  ## Integrate matrix function `f` over a triangle.
  let weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
  let gauss_points = [[1.0/6.0, 1.0/6.0], [2.0/3.0, 1.0/6.0], [1.0/6.0, 2.0/3.0]]
  assert weights.len == gauss_points.len
  let J = jacobian(nodes)
  let detJ = determinant(J)
  for (w, x) in zip(weights, gauss_points):
    result = result + f(x, J) * (w * 0.5) * detJ

proc triSecondOrderQuadratureVec*[N, D](f: proc(x: Vector[D], J: Matrix[D, D]):
    Vector[N], nodes: array[3, array[D, float]]): Vector[N] =
  ## Integrate vector function `f` over a triangle.
  let weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
  let gauss_points = [[1.0/6.0, 1.0/6.0], [2.0/3.0, 1.0/6.0], [1.0/6.0, 2.0/3.0]]
  assert weights.len == gauss_points.len
  let J = jacobian(nodes)
  let detJ = determinant(J)
  for (w, x) in zip(weights, gauss_points):
    result = result + (w * 0.5 * detJ) * f(x, J)

proc tetSecondOrderQuadrature*[N, D](f: proc(x: Vector[D], J: Matrix[D, D]):
    Matrix[N, N], nodes: array[4, array[D, float]]): Matrix[N, N] =
  ## Integrate matrix function `f` over a tetrahedron.
  let weights = [1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0]
  let gauss_points = [
    [0.13819660112501048, 0.13819660112501048, 0.13819660112501048],
    [0.58541019662496852, 0.13819660112501048, 0.13819660112501048],
    [0.13819660112501048, 0.58541019662496852, 0.13819660112501048],
    [0.13819660112501048, 0.13819660112501048, 0.58541019662496852]]
  assert weights.len == gauss_points.len
  let J = jacobian(nodes)
  let detJ = abs(determinant(J))
  for (w, x) in zip(weights, gauss_points):
    result = result + f(x, J) * w * detJ

proc tetSecondOrderQuadratureVec*[N, D](f: proc(x: Vector[D], J: Matrix[D, D]):
    Vector[N], nodes: array[4, array[D, float]]): Vector[N] =
  ## Integrate vector function `f` over a tetrahedron.
  let weights = [1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0]
  let gauss_points = [
    [0.13819660112501048, 0.13819660112501048, 0.13819660112501048],
    [0.58541019662496852, 0.13819660112501048, 0.13819660112501048],
    [0.13819660112501048, 0.58541019662496852, 0.13819660112501048],
    [0.13819660112501048, 0.13819660112501048, 0.58541019662496852]]
  assert weights.len == gauss_points.len
  let J = jacobian(nodes)
  let detJ = abs(determinant(J))
  for (w, x) in zip(weights, gauss_points):
    result = result + (w * detJ) * f(x, J)

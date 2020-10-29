import strutils

type Matrix*[M, N: static[int]] = array[M, array[N, float]]
type Vector*[M: static[int]] = array[M, float]

proc `$`*(matrix: Matrix): string =
  result = "["
  for row in matrix:
    result.add("[")
    for val in row:
      result.add($val & ", ")
    result.removeSuffix({',', ' '})
    result.add("]\n ")
  result.removeSuffix({'\n', ' '})
  result.add("]")

proc rows*(matrix: Matrix): int =
  result = matrix.len

proc cols*(matrix: Matrix): int =
  result = matrix[0].len

proc shape*(matrix: Matrix): (int, int) =
  result = (matrix.rows, matrix.cols)

proc `[]`*(matrix: Matrix, row, col: int): float =
  result = matrix[row][col]

proc `[]=`*(matrix: var Matrix, row, col: int, value: float) =
  matrix[row][col] = value

proc transpose*[M, N](matrix: Matrix[M, N]): Matrix[N, M] =
  ## Compute the transpose of `matrix`.
  for i in 0 ..< result.rows:
    for j in 0 ..< result.cols:
      result[i, j] = matrix[j, i]

proc `+`*(a, b: Matrix): Matrix =
  assert a.cols == b.cols
  assert a.rows == b.rows
  for i in 0 ..< a.rows:
    for j in 0 ..< a.cols:
      result[i, j] = a[i, j] + b[i, j]

proc `+`*[N](a, b: Vector[N]): Vector[N] =
  for i in 0 ..< N:
    result[i] = a[i] + b[i]

proc `-`*[N](a, b: Vector[N]): Vector[N] =
  for i in 0 ..< N:
    result[i] = a[i] - b[i]

proc `*`*[M, N](A: Matrix[M, N], b: Vector[N]): Vector[M] =
  for i in 0 ..< M:
    for j in 0 ..< N:
      result[i] += A[i, j] * b[j]

proc `*`*[M, P, N](a: Matrix[M, P]; b: Matrix[P, N]): Matrix[M, N] =
  for i in 0 ..< result.rows:
    for j in 0 ..< result.cols:
      for k in 0 ..< a.cols:
        result[i, j] = result[i, j] + a[i, k] * b[k, j]
        # would like to do this, but can't get it to compile
        # result[i, j] += a[i, k] * b[k, j]

proc `*`*(a: Matrix, b: float): Matrix =
  for i in 0 ..< a.rows:
    for j in 0 ..< a.cols:
      result[i, j] = b * a[i, j]

proc `*`*[N](a: float, b: Vector[N]): Vector[N] =
  for i in 0 ..< N:
    result[i] = a * b[i]

proc determinant*(a: Matrix[2, 2]): float =
  ## Compute the determinant of a 2x2 matrix `a`.
  result = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]

proc determinant*(a: Matrix[3, 3]): float =
  ## Compute the determinant of a 3x3 matrix `a`.
  result = a[0, 0] * a[1, 1] * a[2, 2] + a[0, 1] * a[1, 2] * a[2, 0] + a[0, 2] *
      a[1, 0] * a[2, 1] - a[0, 2] * a[1, 0] * a[2, 1] - a[0, 1] * a[1, 0] * a[2,
      2] - a[0, 0] * a[1, 2] * a[2, 1]

proc inv*(a: Matrix[2, 2]): Matrix[2, 2] =
  ## Compute the inverse of a 2x2 matrix `a`.
  result = [[a[1, 1], -a[0, 1]],
            [-a[1, 0], a[0, 0]]]
  result = result * (1.0 / determinant(a))

proc inv*(a: Matrix[3, 3]): Matrix[3, 3] =
  ## Compute the inverse of a 3x3 matrix `a`.
  let A = a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]
  let B = a[1, 2] * a[2, 0] - a[1, 0] * a[2, 2]
  let C = a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]

  let D = a[0, 2] * a[2, 1] - a[0, 1] * a[2, 2]
  let E = a[0, 0] * a[2, 2] - a[0, 2] * a[2, 0]
  let F = a[0, 1] * a[2, 0] - a[0, 0] * a[2, 1]

  let G = a[0, 1] * a[1, 2] - a[0, 2] * a[1, 1]
  let H = a[0, 2] * a[1, 0] - a[0, 0] * a[1, 2]
  let I = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]

  result = [[A, D, G],
            [B, E, H],
            [C, F, I]] * (1.0 / determinant(a))

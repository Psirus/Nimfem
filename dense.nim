import math
import strutils

type Vector* = seq[float] ## A dynamic vector
type Matrix* = seq[seq[float]] ## A dynamic matrix

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

proc newMatrix*(rows, cols: int): Matrix =
  result = newSeq[seq[float]](rows)
  for i in 0 ..< rows:
    result[i] = newSeq[float](cols)

proc newMatrix*(shape: (int, int)): Matrix =
  result = newMatrix(shape[0], shape[1])

proc newVector*(size: int): Vector =
  result = newSeq[float](size)

proc size*(vector: Vector): int =
  result = vector.len

proc shape*(matrix: Matrix): (int, int) =
  result = (matrix.rows, matrix.cols)

proc `[]`*(matrix: Matrix, row, col: int): float =
  result = matrix[row][col]

proc `[]=`*(matrix: var Matrix, row, col: int, value: float) =
  matrix[row][col] = value

proc transpose*(matrix: Matrix): Matrix =
  result = newMatrix(matrix.cols, matrix.rows)
  for i in 0 ..< result.rows:
    for j in 0 ..< result.cols:
      result[i, j] = matrix[j, i]

proc `+`*(a, b: Matrix): Matrix =
  assert a.cols == b.cols
  assert a.rows == b.rows
  result = newMatrix(a.rows, a.cols)
  for i in 0 ..< a.rows:
    for j in 0 ..< a.cols:
      result[i, j] = a[i, j] + b[i, j]

proc `+`*(a, b: Vector): Vector =
  assert a.size == b.size
  result = newVector(a.size)
  for i in 0 ..< a.size:
    result[i] = a[i] + b[i]

proc `-`*(a, b: Vector): Vector =
  assert a.size == b.size
  result = newVector(a.size)
  for i in 0 ..< a.size:
    result[i] = a[i] - b[i]

proc `*`*(A: Matrix, b: Vector): Vector =
  assert A.cols == b.size
  result = newVector(b.size)
  for i in 0 ..< result.size:
    for j in 0 ..< A.cols:
      result[i] += A[i, j] * b[j]

proc `*`*(a, b: Matrix): Matrix =
  assert a.cols == b.rows
  result = newMatrix(a.rows, b.cols)
  for i in 0 ..< result.rows:
    for j in 0 ..< result.cols:
      for k in 0 ..< a.cols:
        result[i, j] = result[i, j] + a[i, k] * b[k, j]
        # would like to do this, but can't get it to compile
        # result[i, j] += a[i, k] * b[k, j]

proc `*`*(a: Matrix, b: float): Matrix =
  result = newMatrix(a.shape)
  for i in 0 ..< a.rows:
    for j in 0 ..< a.cols:
      result[i, j] = b * a[i, j]

proc `*`*(a: float, b: Vector): Vector =
  result = newVector(b.size)
  for i in 0 ..< b.size:
      result[i] = a * b[i]

proc determinant*(a: Matrix): float =
  # currently only for 2x2
  assert a.shape == (2, 2)
  result = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]

proc inv*(a: Matrix): Matrix =
  # currently only for 2x2
  assert a.shape == (2, 2)
  result = @[@[ a[1, 1], -a[0, 1]],
             @[-a[1, 0],  a[0, 0]]]
  result = result * (1.0/determinant(a))

proc dot*(a, b: Vector): float =
  assert a.size == b.size
  for i in 0..<a.size:
    result += a[i] * b[i]

proc norm*(a: Vector): float =
  result = sqrt(dot(a, a))

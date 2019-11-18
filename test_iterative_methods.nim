import sequtils
import fenv

import sparse
import iterative_methods

let eps = 2 * epsilon float

proc allClose(a, b: DynamicVector): bool =
  var tmp = newVector(a.size)
  for i in 0..<a.size:
    tmp[i] = abs(a[i] - b[i])
  result = all(tmp, proc (x: float): bool = return x < eps)

block:
  var A: SparseMatrix
  A.ia = @[0, 2, 4]
  A.ja = @[0, 1, 0, 1]
  A.aa = @[4.0, 1.0, 1.0, 3.0]

  let b = @[1.0, 2.0]
  let u = conjugate_gradient(A, b)
  doAssert allClose(u, @[1.0/11.0, 7.0/11.0])

block:
  var A: SparseMatrix
  A.ia = @[0, 1, 2, 3, 4]
  A.ja = @[0, 1, 2, 3]
  A.aa = @[1.0, 2.0, 3.0, 4.0]

  var b = @[1.0, 1.0, 1.0, 1.0]
  let u = conjugate_gradient(A, b)
  doAssert allClose(u, @[1.0, 0.5, 1.0/3.0, 0.25])

import sequtils
import fenv

import ../src/sparse
import ../src/iterative_methods

let eps = 2 * epsilon float

proc allClose(a, b: DynamicVector, atol = eps, rtol = eps): bool =
  var tmp = newSeq[bool](a.size)
  for i in 0..<a.size:
    tmp[i] = abs(a[i] - b[i]) <= (atol + rtol * abs(b[i]))
  result = all(tmp, proc (x: bool): bool = return x)

block:
  var A: SparseMatrix
  A.ia = @[0, 2, 4]
  A.ja = @[0, 1, 0, 1]
  A.aa = @[4.0, 1.0, 1.0, 3.0]

  let b = @[1.0, 2.0]
  let u = conjugate_gradient(A, b)
  doAssert allClose(u, @[1.0/11.0, 7.0/11.0])

block:
  let Ai = @[0, 1, 2, 3]
  let Aj = @[0, 1, 2, 3]
  let Ax = @[1.0, 2.0, 3.0, 4.0]
  let A = toCSR(Ai, Aj, Ax)

  var b = @[1.0, 1.0, 1.0, 1.0]
  let u = conjugate_gradient(A, b)
  doAssert allClose(u, @[1.0, 0.5, 1.0/3.0, 0.25])

block:
  let Ai = @[0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4]
  let Aj = @[0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4]
  let Ax = @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
  let A = toCSR(Ai, Aj, Ax)

  let ILUx = @[1.0, 2.0, 3.0, 4.0, -1.0, 6.0, 7.0, -4.0, 9.0, 1.428571428571429, 1.671428571428572e+01, 12.0]
  let reference_ilu = toCSR(Ai, Aj, ILUx)

  let computed_ilu = incomplete_lu(A)

  doAssert reference_ilu.ia == computed_ilu.ia
  doAssert reference_ilu.ja == computed_ilu.ja
  doAssert allClose(reference_ilu.aa, computed_ilu.aa)

block:
  let Ai = @[0, 1, 2, 3, 0, 1, 2, 2]
  let Aj = @[0, 1, 2, 3, 2, 0, 0, 1]
  let Ax = @[1.0, 2.0, 3.0, 4.0, 4.0, 1.0, 1.0, 2.0]
  let A = toCSR(Ai, Aj, Ax)

  let ILUx = @[1.0, 2.0, -1.0, 4.0, 4.0, 1.0, 1.0, 1.0]
  let reference_ilu = toCSR(Ai, Aj, ILUx)

  let computed_ilu = incomplete_lu(A)

  doAssert reference_ilu.ia == computed_ilu.ia
  doAssert reference_ilu.ja == computed_ilu.ja
  doAssert allClose(reference_ilu.aa, computed_ilu.aa)

block:
  let Ai = @[0, 0, 0, 1, 1, 1, 2, 2, 2]
  let Aj = @[0, 1, 2, 0, 1, 2, 0, 1, 2]
  let Ax = @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 2.0]
  let A = toCSR(Ai, Aj, Ax)

  let ILUx = @[1.0, 2.0, 3.0, 4.0, -3.0, -6.0, 7.0, 2.0, -7.0]
  let reference_ilu = toCSR(Ai, Aj, ILUx)

  let computed_ilu = incomplete_lu(A)

  doAssert reference_ilu.ia == computed_ilu.ia
  doAssert reference_ilu.ja == computed_ilu.ja
  doAssert allClose(reference_ilu.aa, computed_ilu.aa)

block:
  let Ai = @[0, 1, 2, 3, 0, 1, 2, 2]
  let Aj = @[0, 1, 2, 3, 2, 0, 0, 1]
  let Ax = @[1.0, 2.0, 3.0, 4.0, 4.0, 1.0, 1.0, 2.0]
  let A = toCSR(Ai, Aj, Ax)

  let b = @[1.0, 2.0, 3.0, 4.0]
  let P = incomplete_lu(A)
  let x = solve_ilu(P, b)

  doAssert allClose(x, @[5.0, 0.5, -1.0, 1.0])

block:
  var A: SparseMatrix
  A.ia = @[0, 2, 4]
  A.ja = @[0, 1, 0, 1]
  A.aa = @[4.0, 1.0, 1.0, 3.0]

  let b = @[1.0, 2.0]
  let P = incomplete_lu(A)
  let u = preconditioned_cg(A, P, b)
  doAssert allClose(u, @[1.0/11.0, 7.0/11.0])

block:
  let Ai = @[0, 1, 2, 3]
  let Aj = @[0, 1, 2, 3]
  let Ax = @[1.0, 2.0, 3.0, 4.0]
  let A = toCSR(Ai, Aj, Ax)

  var b = @[1.0, 1.0, 1.0, 1.0]
  let P = incomplete_lu(A)
  let u = preconditioned_cg(A, P, b)
  doAssert allClose(u, @[1.0, 0.5, 1.0/3.0, 0.25])

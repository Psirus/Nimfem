import sequtils
import fenv

import ../src/dense
import ../src/tetrahedron

let eps = 2 * epsilon float

proc allClose(a, b: Vector, atol = eps, rtol = eps): bool =
  var tmp = newSeq[bool](a.len)
  for i in 0..<a.len:
    tmp[i] = abs(a[i] - b[i]) <= (atol + rtol * abs(b[i]))
  result = all(tmp, proc (x: bool): bool = return x)

block:
  let x = [0.0, 0.0, 0.0]
  let N = shapeFunctions(x)
  doAssert allClose(N, [1.0, 0.0, 0.0, 0.0])

block:
  let nodes = [[0.0, 0.0, 0.0],
               [2.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [0.0, 0.0, 1.5]]
  let J = jacobian(nodes)
  doAssert determinant(J) == 3.0

block:
  let nodes = [[0.0, 0.0, 0.0],
               [1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [0.0, 0.0, 1.0]]
  let invJ = inv(jacobian(nodes))
  let B = kelvinStrain(invJ)

  var u = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  var strain = B*u
  doAssert strain == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

  u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
  strain = B*u
  doAssert strain == [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

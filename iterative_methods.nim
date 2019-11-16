import dense
import sparse

proc conjugate_gradient*(A: SparseMatrix, b: Vector): Vector =
  result = b
  var R = b - A * result
  var P = R
  var m = 0 
  let tol = 1e-9
  while norm(R) > tol:
    let Ap = A * P
    let alpha = dot(R, R) / dot(P, Ap)
    result = result + alpha * P
    let scale = dot(R, R)
    R = R - alpha * Ap
    let beta = dot(R, R) / scale
    P = R + beta * P
    m += 1
    echo "res = ", norm(R)

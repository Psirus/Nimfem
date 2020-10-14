import math
import sparse

proc conjugate_gradient*(A: SparseMatrix, b: DynamicVector): DynamicVector =
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
    # echo "res = ", norm(R)
  # echo m

proc isNan(x: float): bool =
  result = classify(x) == fcNan

# Saad
proc incomplete_lu*(A: SparseMatrix): SparseMatrix =
  result = A
  for k in 0 ..< A.rows-1:
    let d = 1.0 / result.getEntry(k,k)
    for i in (k+1) ..< A.rows:
      let entry = result.getEntry(i, k)
      if not isNan(entry):
        let e = d * entry
        result.setEntry(i, k, e)
        for j in (k+1) ..< A.rows:
          let A_ij = result.getEntry(i, j)
          let A_kj = result.getEntry(k, j)
          if (not isNan(A_ij)) and (not isNan(A_kj)):
            result.setEntry(i, j, A_ij - e * A_kj)

proc solve_ilu*(P: SparseMatrix, b: DynamicVector): DynamicVector =
  assert b.len == P.cols, "Preconditioner and RHS vector don't match in size."
  result = b
  # forward substitution 
  for i in 0 ..< P.rows:
    var sum = 0.0
    let nz_bounds = P.nonzero_bounds_row(i)
    for j_idx in nz_bounds[0] .. nz_bounds[1]:
      let j = P.ja[j_idx]
      if j >= i:
        break
      sum += P.aa[j_idx] * result[j]
    result[i] -= sum

  # backward substitution
  for i in countdown(P.rows-1, 0):
    var sum = 0.0
    let nz_bounds = P.nonzero_bounds_row(i)
    for j_idx in nz_bounds[0] .. nz_bounds[1]:
      let j = P.ja[j_idx]
      if j <= i:
        continue
      sum += P.aa[j_idx] * result[j]
    result[i] = (result[i] - sum) / P.getEntry(i, i)

proc preconditioned_cg*(A: SparseMatrix, C: SparseMatrix, b: DynamicVector): DynamicVector =
  result = b
  var R = b - A * result
  var z = solve_ilu(C, R)
  var P = z
  var m = 0
  let tol = 1e-9
  while norm(R) > tol:
    let Ap = A*P
    let zR = dot(z, R)
    let alpha = zR / dot(Ap, P)
    result = result + alpha * P
    R = R - alpha * Ap
    z = solve_ilu(C, R)
    let beta = dot(z, R) / zR
    P = z + beta * P
    m += 1
    # echo "res = ", norm(R)
  # echo m

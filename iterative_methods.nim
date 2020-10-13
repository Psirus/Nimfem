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
    echo "res = ", norm(R)
  echo m

proc isNan(x: float): bool =
  result = classify(x) == fcNan

# Ern & Guermond
proc incomplete_lu*(A: SparseMatrix): SparseMatrix =
  result = A
  for i in 0 ..< A.rows:
    for k in 0 ..< i:
      let entry = A.getEntry(i, k)
      if not isNan(entry):
        let e = entry / A.getEntry(k, k)
        result.setEntry(i, k, e)
        for j in k+1 .. A.rows:
          let A_ij = A.getEntry(i, j)
          let A_kj = A.getEntry(k, j)
          if (not isNan(A_ij)) and (not isNan(A_kj)):
            let A_ik = A.getEntry(i, k)
            result.setEntry(i, j, A_ij - A_ik * A_kj)

# CFD Online
proc incomplete_lu2*(A: SparseMatrix): SparseMatrix =
  result = A
  for r in 0 ..< A.rows-1:
    let d = 1.0 / A.getEntry(r, r)
    for i in (r+1) ..< A.rows:
      let entry = A.getEntry(i, r)
      if not isNan(entry):
        let e = d * entry
        result.setEntry(i, r, e)
        for j in r+1 ..< A.rows:
          let A_ij = A.getEntry(i, j)
          let A_rj = A.getEntry(r, j)
          if (not isNan(A_ij)) and (not isNan(A_rj)):
            result.setEntry(i, j, A_ij - e * A_rj)

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
    echo "res = ", norm(R)

  echo m

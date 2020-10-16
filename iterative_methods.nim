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

# Saad, Algorithm 10.3
proc incomplete_lu*(A: SparseMatrix): SparseMatrix =
  result = A
  for i in 1 ..< A.rows:
    let nz_bounds = result.nonzero_bounds_row(i)
    for k_idx in nz_bounds[0] .. nz_bounds[1]:
      let k = result.ja[k_idx]
      if k >= i:
        break
      let e = result.aa[k_idx] / result.getEntry(k, k)
      result.aa[k_idx] = e
      for j_idx in nz_bounds[0] .. nz_bounds[1]:
        let j = result.ja[j_idx]
        if j < (k+1):
          continue
        let kj_idx = result.getIndex(k, j)
        if kj_idx != -1:
          result.aa[j_idx] -= e * result.aa[kj_idx]

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

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

proc incomplete_lu*(A: SparseMatrix): SparseMatrix =
  result = A
  for i in 1 ..< A.rows:
    for k in 0 ..< i:
      let entry = A.getEntry(i, k)
      if bool(entry):
        let e = entry / A.getEntry(k, k)
        result.setEntry(i, k, e)
        for j in k+1 .. A.rows:
          let A_ij = A.getEntry(i, j)
          let A_kj = A.getEntry(k, j)
          if bool(A_ij) and bool(A_kj):
            let A_ik = A.getEntry(i, k)
            result.setEntry(i, j, A_ij - A_ik * A_kj)

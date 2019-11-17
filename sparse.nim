import algorithm
import sequtils
import fenv
import dense
let eps = 2 * epsilon float

type
  SparseMatrix* = object
    ia*: seq[int]
    ja*: seq[int]
    aa*: seq[float]

proc setDiagonalRows*(A: var SparseMatrix, rows: seq[int]) =
  for row in rows:
    let nnz_in_row = A.ia[row+1] - A.ia[row]
    let lower = A.ia[row]
    let upper = A.ia[row] + nnz_in_row - 1
    for idx in row+1..<A.ia.len:
      A.ia[idx] -= nnz_in_row - 1
    A.ja.delete(lower, upper)
    A.aa.delete(lower, upper)

    A.ja.insert(row, lower)
    A.aa.insert(1.0, lower)

proc `*`*(A: SparseMatrix, b: Vector): Vector =
  result = newVector(b.size)
  for i in 0..<A.ia.len-1:
    for idx in A.ia[i]..<A.ia[i+1]:
      let col = A.ja[idx]
      result[i] += A.aa[idx] * b[col]

proc toDense*(A: SparseMatrix): Matrix =
  let rows = A.ia.len - 1
  let cols = max(A.ja) + 1
  result = newMatrix(rows, cols)
  for i in 0..<rows:
    for idx in A.ia[i]..<A.ia[i+1]:
      let j = A.ja[idx]
      result[i, j] = A.aa[idx]
    
proc `$`*(A: SparseMatrix): string =
  let A_dense = toDense(A)
  result = $A_dense

proc removeZeros(A: var SparseMatrix) =
  var row_end = 0
  var nnz = 0
  for i in 0..<A.ia.high:
    var jj = row_end
    row_end = A.ia[i+1]
    while jj < row_end:
      let j = A.ja[jj]
      let val = A.aa[jj]
      if abs(val) > eps:
        A.ja[nnz] = j
        A.aa[nnz] = val
        nnz += 1
      jj += 1
    A.ia[i+1] = nnz

  A.ja.delete(nnz, A.ja.high)
  A.aa.delete(nnz, A.aa.high)

proc sortCols(A: var SparseMatrix) =
  for i in 0..<A.ia.high:
    let row_start = A.ia[i]
    let row_end = A.ia[i+1]
    var temp: seq[tuple[col: int, val: float]]
    for idx in row_start..<row_end:
      temp.add((A.ja[idx], A.aa[idx]))
    temp = temp.sortedByIt(it.col)

    var n = 0
    for idx in row_start..<row_end:
      A.ja[idx] = temp[n][0]
      A.aa[idx] = temp[n][1]
      n += 1

proc sumDuplicates(A: var SparseMatrix) =
  var row_end = 0
  var nnz = 0
  for i in 0..<A.ia.high:
    var jj = row_end
    row_end = A.ia[i+1]
    while jj < row_end:
      let j = A.ja[jj]
      var val = A.aa[jj]
      jj += 1
      while (jj < row_end) and (A.ja[jj] == j):
        val += A.aa[jj]
        jj += 1
      A.ja[nnz] = j
      A.aa[nnz] = val
      nnz += 1
    A.ia[i+1] = nnz

  A.ja.delete(nnz, A.ja.high)
  A.aa.delete(nnz, A.aa.high)

proc cumSum(xs: var seq[int]) =
  for i in 1..<xs.len:
    xs[i] += xs[i-1]

proc toCSR*(Ai, Aj: seq[int], Ax: seq[float]): SparseMatrix =
  let n_rows = max(Ai) + 1
  result.ia = newSeq[int](n_rows + 1)
  for row in Ai:
    result.ia[row+1] += 1
  cumSum(result.ia)

  result.ja = newSeq[int](Aj.len)
  result.aa = newSeq[float](Ax.len)
  for i in 0..<Aj.len:
    let row = Ai[i]
    let dest = result.ia[row]
    result.ja[dest] = Aj[i]
    result.aa[dest] = Ax[i]
    result.ia[row] += 1

  var last = 0
  for i in 0..n_rows:
    let tmp = result.ia[i]
    result.ia[i] = last
    last = tmp

  removeZeros(result)
  sortCols(result)
  sumDuplicates(result)

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

proc fromTripletList*(triplets: seq[(int, int, float)]): SparseMatrix =
  var nrows = 0
  var nnz_per_row: seq[int]
  var col_val_pairs: seq[seq[tuple[col: int, val: float]]]
  # remove entries close to zero, add up duplicate entries
  var indices: seq[(int, int)]
  var values: seq[float]
  for (row, col, val) in triplets:
    if abs(val) < eps:
      continue
    let idx = indices.find((row, col))
    if idx == -1:
      indices.add((row, col))
      values.add(val)
    else:
      values[idx] += val

  for (indices, val) in zip(indices, values):
    # this if block is just dynamic growth, not relevant for the algorithm
    # therefore we could be faster if we knew the number of rows
    let (row, col) = indices
    if row >= nrows:
      let growth = row - nrows + 1
      nnz_per_row.add(newSeq[int](growth))
      col_val_pairs.add(newSeq[seq[(int, float)]](growth))
      nrows = row + 1
    # count number of nnz in each row
    nnz_per_row[row] += 1
    col_val_pairs[row].add((col, val))

  # construct the ia array
  result.ia = newSeq[int](nnz_per_row.len + 1)
  result.ia[0] = 1
  for i, nnz in nnz_per_row.pairs:
    result.ia[i + 1] = result.ia[i] + nnz

  # construct the ja, aa arrays
  result.ja = newSeqOfCap[int](result.ia.high - 1)
  result.aa = newSeqOfCap[float](result.ia.high - 1)
  for row in col_val_pairs:
    let col_val_pairs  = row.sortedByIt(it.col)
    for (col, val) in col_val_pairs:
      result.ja.add(col)
      result.aa.add(val)


proc setDiagonalRows*(A: var SparseMatrix, rows: seq[int]) =
  for row in rows:
    let nnz_in_row = A.ia[row+1] - A.ia[row]
    let lower = A.ia[row] - 1
    let upper = A.ia[row] + nnz_in_row - 2
    for idx in row+1..<A.ia.len:
      A.ia[idx] -= nnz_in_row - 1
    A.ja.delete(lower, upper)
    A.aa.delete(lower, upper)

    A.ja.insert(row, lower)
    A.aa.insert(1.0, lower)

proc `*`*(A: SparseMatrix, b: Vector): Vector =
  result = newVector(b.size)
  for i in 0..<A.ia.len-1:
    for idx in A.ia[i]-1..<A.ia[i+1]-1:
      let col = A.ja[idx]
      result[i] += A.aa[idx] * b[col]

proc toDense*(A: SparseMatrix): Matrix =
  let rows = A.ia.len - 1
  result = newMatrix(rows, rows)
  for i in 0..<rows:
    for idx in A.ia[i]-1..<A.ia[i+1]-1:
      let j = A.ja[idx]
      result[i, j] = A.aa[idx]
    
proc `$`*(A: SparseMatrix): string =
  let A_dense = toDense(A)
  result = $A_dense

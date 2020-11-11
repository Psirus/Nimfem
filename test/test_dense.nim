import sequtils
import fenv
import ../src/dense

let eps = 2 * epsilon float

proc allClose(a, b: Matrix, atol = eps, rtol = eps): bool =
  if (a.shape != b.shape):
    return false
  var tmp = newSeq[bool](a.rows * a.cols)
  for i in 0..<a.rows:
    for j in 0..<a.cols:
      tmp[a.rows*i + j] = abs(a[i,j] - b[i,j]) <= (atol + rtol * abs(b[i, j]))
  result = all(tmp, proc (x: bool): bool = return x)

let a = [[1.0, 1.0, 1.0, 1.0],
          [2.0, 4.0, 8.0, 16.0],
          [3.0, 9.0, 27.0, 81.0],
          [4.0, 16.0, 64.0, 256.0]]

let b = [[4.0, -3.0, 4/3.0, -1/4.0],
          [-13/3.0, 19/4.0, -7/3.0, 11/24.0],
          [3/2.0, -2.0, 7/6.0, -1/4.0],
          [ -1/6.0, 1/4.0, -1/6.0, 1/24.0]]

echo dense.`$`(a)
echo dense.`$`(b)
echo dense.`$`(a * b)
echo dense.`$`(b * a)

block:
  # just np.random.rand(3,3) and np.linalg.inv results
  let A = [[0.6317672566505066, 0.6468079262002876, 0.7189160871021079],
           [0.2926953841937896, 0.4788163688370418, 0.2394338892261617],
           [0.2148476885155944, 0.927548744504706 , 0.1671958197892681]]
  let expected = [[ -4.2895466873950605,  16.873207611279245 , -5.719004239163235 ],
                  [  0.0756360363273769,  -1.4747016490599907,   1.7866331130693196],
                  [  5.092490237524678 , -13.501007338984667 ,   3.418300416825433 ]]
  doAssert allClose(inv(A), expected)

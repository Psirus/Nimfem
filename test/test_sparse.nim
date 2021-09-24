import unittest

import ../src/sparse

suite "Sparse Matrices":
  test "To CSR":
    let Ai = @[0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4]
    let Aj = @[0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4]
    let Ax = @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    var m = toCSR(Ai, Aj, Ax)
    check m.ia == @[0, 2, 5, 9, 11, 12]
    check m.ja == @[0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4]
    check m.aa == @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]

    m.setDiagonalRows(@[2])
    check m.ia == @[0, 2, 5, 6, 8, 9]
    check m.ja == @[0, 3, 0, 1, 3, 2, 2, 3, 4]
    check m.aa == @[1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 10.0, 11.0, 12.0]

    var b = m * @[1.0, 2.0, 3.0, 4.0, 5.0]
    check b == @[9.0, 31.0, 3.0, 74.0, 60.0]

  test "Sum duplicates":
    let Ai = @[0, 2, 1, 0, 2, 1, 0, 2]
    let Aj = @[3, 5, 1, 2, 7, 4, 7, 5]
    let Ax = @[3.0, 7.0, 2.0, 1.0, 1.0, 1.0, 5.0, 3.0]
    var m = toCSR(Ai, Aj, Ax)

    check m.ia == @[0, 3, 5, 7]
    check m.ja == @[2, 3, 7, 1, 4, 5, 7]
    check m.aa == @[1.0, 3.0, 5.0, 10.0, 1.0]

    skip()

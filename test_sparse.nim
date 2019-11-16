import sparse

block:
  let triplets = @[(0, 0, 1.0), (0, 3, 2.0), (1, 0, 3.0), (1, 1, 4.0),
                   (1, 3, 5.0), (2, 0, 6.0), (2, 2, 7.0), (2, 3, 8.0),
                   (2, 4, 9.0), (3, 2, 10.0), (3, 3, 11.0), (4, 4, 12.0)]
  var m = fromTripletList(triplets)
  doAssert m.ia == @[1, 3, 6, 10, 12, 13]
  doAssert m.ja == @[0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4]
  doAssert m.aa == @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]

  m.setDiagonalRows(@[2])
  doAssert m.ia == @[1, 3, 6, 7, 9, 10]
  doAssert m.ja == @[0, 3, 0, 1, 3, 2, 2, 3, 4]
  doAssert m.aa == @[1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 10.0, 11.0, 12.0]

block:
  let triplets = @[(0, 3, 3.0), (2, 5, 7.0), (0, 2, 1.0), (2, 7, 1.0), (0, 7, 5.0)]
  let m = fromTripletList(triplets)
  doAssert m.ia == @[1, 4, 4, 6]
  doAssert m.ja == @[2, 3, 7, 5, 7]
  doAssert m.aa == @[1.0, 3.0, 5.0, 7.0, 1.0]

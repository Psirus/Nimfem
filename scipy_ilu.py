from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spilu

A = csc_matrix([[1.0, 0.0, 4.0, 0.0],
 [1.0, 2.0, 0.0, 0.0],
 [1.0, 2.0, 3.0, 0.0],
 [0.0, 0.0, 0.0, 4.0]])

B = spilu(A)
from IPython import embed; embed()

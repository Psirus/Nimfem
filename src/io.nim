import strformat

import mesh
import sparse

proc writeVTK*[Element](filename: string, mesh: Mesh[Element], u: DynamicVector) =
  ## Write VTK output to `filename`.
  var vtkFile = open(filename, FileMode.fmWrite)
  # Header
  vtkFile.write("# vtk DataFile Version 1.0\n")
  vtkFile.write("Nimfem output\n")
  vtkFile.write("ASCII\n\n")

  # nodes
  let num_nodes = mesh.nodes.len
  vtkFile.write("DATASET UNSTRUCTURED_GRID\n")
  vtkFile.write("POINTS {num_nodes} float\n".fmt)
  for node in mesh.nodes:
    when Element.dim == 2:
      vtkFile.write("{node[0]} {node[1]} 0.0\n".fmt)
    elif Element.dim == 3:
      vtkFile.write("{node[0]} {node[1]} {node[2]}\n".fmt)
  vtkFile.write("\n")

  # cells
  let num_cells = mesh.connectivity.len
  let nodes_per_cell = Element.num_nodes
  vtkFile.write("CELLS {num_cells} {(nodes_per_cell+1)*num_cells}\n".fmt)
  for cell in mesh.connectivity:
    when Element.num_nodes == 3:
      vtkFile.write("3 {cell[0]} {cell[1]} {cell[2]}\n".fmt)
    elif Element.num_nodes == 4:
      vtkFile.write("4 {cell[0]} {cell[1]} {cell[2]} {cell[3]}\n".fmt)
  vtkFile.write("\n")
  vtkFile.write("CELL_TYPES {num_cells}\n".fmt)
  for _ in 0..<num_cells:
    when Element.shape == triangle:
      vtkFile.write("5\n")
    elif Element.shape == tetrahedron:
      vtkFile.write("10\n")
  vtkFile.write("\n")

  ## scalar data
  #vtkFile.write("POINT_DATA {num_nodes}\n".fmt)
  #vtkFile.write("SCALARS u float\n")
  #vtkFile.write("LOOKUP_TABLE default\n")
  #for val in u:
  #  vtkFile.write("{val}\n".fmt)

  # vector data
  vtkFile.write("POINT_DATA {num_nodes}\n".fmt)
  vtkFile.write("VECTORS u float\n")
  for i in countup(0, u.len-1, 3):
    vtkFile.write("{u[i]} {u[i+1]} {u[i+2]}\n".fmt)

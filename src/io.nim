import strformat

import mesh
import sparse

proc writeVTK*(filename: string, mesh: Mesh, u: DynamicVector) =
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
    vtkFile.write("{node[0]} {node[1]} 0.0\n".fmt)
  vtkFile.write("\n")

  # cells
  let num_cells = mesh.connectivity.len
  vtkFile.write("CELLS {num_cells} {4*num_cells}\n".fmt)
  for cell in mesh.connectivity:
    vtkFile.write("3 {cell[0]} {cell[1]} {cell[2]}\n".fmt)
  vtkFile.write("\n")
  vtkFile.write("CELL_TYPES {num_cells}\n".fmt)
  for _ in 0..<num_cells:
    vtkFile.write("5\n")

  # data
  vtkFile.write("POINT_DATA {num_nodes}\n".fmt)
  vtkFile.write("SCALARS u float\n")
  vtkFile.write("LOOKUP_TABLE default\n")
  for val in u:
    vtkFile.write("{val}\n".fmt)

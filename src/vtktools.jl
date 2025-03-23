# Create and return VTK grids on enriched uniform nodes
# that are ready to be filled with reinterpolated data (vtu version)
function build_vtk_grids(::Val{:vtu}, mesh::TreeMesh, nodes, n_visnodes, verbose,
                         output_directory, is_datafile, filename, reinterpolate::Val{true})
  coordinates, levels, center_level_0, length_level_0 = extract_mesh_information(mesh)

  # Extract number of spatial dimensions
  ndims_ = size(coordinates, 1)

  # Calculate VTK points and cells
  verbose && println("| Preparing VTK cells...")
  if is_datafile
    @timeit "prepare VTK cells (node data)" begin
      vtk_points, vtk_cells = calc_vtk_points_cells(Val(ndims_), coordinates, levels,
                                                    center_level_0, length_level_0, n_visnodes)
    end
  end

  # Prepare VTK points and cells for celldata file
  @timeit "prepare VTK cells (cell data)" begin
    vtk_celldata_points, vtk_celldata_cells = calc_vtk_points_cells(Val(ndims_), coordinates,
                                                                    levels, center_level_0,
                                                                    length_level_0, 1)
  end

  # Determine output file names
  base, _ = splitext(splitdir(filename)[2])
  vtk_filename = joinpath(output_directory, base)
  vtk_celldata_filename = vtk_filename * "_celldata"

  # Open VTK files
  verbose && println("| Building VTK grid...")
  if is_datafile
    @timeit "build VTK grid (node data)" vtk_nodedata = vtk_grid(vtk_filename, vtk_points,
                                                                 vtk_cells)
  else
    vtk_nodedata = nothing
  end
  @timeit "build VTK grid (cell data)" vtk_celldata = vtk_grid(vtk_celldata_filename,
                                                                vtk_celldata_points,
                                                                vtk_celldata_cells)

  return vtk_nodedata, vtk_celldata
end


# Create and return VTK grids with nodes on which the raw was created
# ready to be filled with raw data (vtu version)
function build_vtk_grids(::Val{:vtu}, mesh::TreeMesh,
                         nodes, n_visnodes, verbose,
                         output_directory, is_datafile, filename, reinterpolate::Val{false})

  @timeit "prepare coordinate information" node_coordinates = calc_node_coordinates(mesh, nodes, n_visnodes)

  # Calculate VTK points and cells
  verbose && println("| Preparing VTK cells...")
  if is_datafile
    @timeit "prepare VTK cells (node data)" begin
      vtk_points, vtk_cells = calc_vtk_points_cells(node_coordinates)
    end
  end

  # Prepare VTK points and cells for celldata file
  @timeit "prepare VTK cells (cell data)" begin
    vtk_celldata_points, vtk_celldata_cells = calc_vtk_points_cells(node_coordinates)
  end

  # Determine output file names
  base, _ = splitext(splitdir(filename)[2])
  vtk_filename = joinpath(output_directory, base)
  vtk_celldata_filename = vtk_filename * "_celldata"

  # Open VTK files
  verbose && println("| Building VTK grid...")
  if is_datafile
    @timeit "build VTK grid (node data)" vtk_nodedata = vtk_grid(vtk_filename, vtk_points,
                                                                 vtk_cells)
  else
    vtk_nodedata = nothing
  end
  @timeit "build VTK grid (cell data)" vtk_celldata = vtk_grid(vtk_celldata_filename,
                                                                vtk_celldata_points,
                                                                vtk_celldata_cells)

  return vtk_nodedata, vtk_celldata
end


# Create and return VTK grids on enriched uniform nodes
# that are ready to be filled with reinterpolated data (vti version)
function build_vtk_grids(::Val{:vti}, mesh::TreeMesh, nodes, n_visnodes, verbose,
                         output_directory, is_datafile, filename, reinterpolate::Val{true})

  coordinates, levels, center_level_0, length_level_0 = extract_mesh_information(mesh)

  # Extract number of spatial dimensions
  ndims_ = size(coordinates, 1)

  # Prepare VTK points and cells for celldata file
  @timeit "prepare VTK cells (cell data)" begin
    vtk_celldata_points, vtk_celldata_cells = calc_vtk_points_cells(Val(ndims_), coordinates,
                                                                    levels, center_level_0,
                                                                    length_level_0, 1)
  end

  # Determine output file names
  base, _ = splitext(splitdir(filename)[2])
  vtk_filename = joinpath(output_directory, base)
  vtk_celldata_filename = vtk_filename * "_celldata"

  # Open VTK files
  verbose && println("| Building VTK grid...")
  if is_datafile
    # Determine level-wise resolution
    max_level = maximum(levels)
    resolution = n_visnodes * 2^max_level
    origin = @. center_level_0 - 0.5 * length_level_0

    Nx = Ny = resolution + 1
    dx = dy = length_level_0/resolution
    spacing = (dx, dy)
    @timeit "build VTK grid (node data)" vtk_nodedata = vtk_grid(vtk_filename, Nx, Ny,
                                                                 origin=tuple(origin...),
                                                                 spacing=spacing)
  else
    vtk_nodedata = nothing
  end
  @timeit "build VTK grid (cell data)" vtk_celldata = vtk_grid(vtk_celldata_filename,
                                                      vtk_celldata_points,
                                                      vtk_celldata_cells)

  return vtk_nodedata, vtk_celldata
end


# Create and return VTK grids that are ready to be filled with data
# (StructuredMesh/UnstructuredMesh2D/P4estMesh version).
# Routine is agnostic with respect to reinterpolation.
function build_vtk_grids(::Val{:vtu},
                         mesh::Union{StructuredMesh, UnstructuredMesh2D, P4estMesh, T8codeMesh},
                         nodes, n_visnodes, verbose, output_directory, is_datafile, filename,
                         reinterpolate::Union{Val{true}, Val{false}})

  @timeit "prepare coordinate information" node_coordinates = calc_node_coordinates(mesh, nodes, n_visnodes)

  # Calculate VTK points and cells
  verbose && println("| Preparing VTK cells...")
  if is_datafile
    @timeit "prepare VTK cells (node data)" begin
      vtk_points, vtk_cells = calc_vtk_points_cells(node_coordinates)
    end
  end

  # Prepare VTK points and cells for celldata file
  @timeit "prepare VTK cells (cell data)" begin
    vtk_celldata_points, vtk_celldata_cells = calc_vtk_points_cells(node_coordinates)
  end

  # Determine output file names
  base, _ = splitext(splitdir(filename)[2])
  vtk_filename = joinpath(output_directory, base)
  vtk_celldata_filename = vtk_filename * "_celldata"

  # Open VTK files
  verbose && println("| Building VTK grid...")
  if is_datafile
    @timeit "build VTK grid (node data)" vtk_nodedata = vtk_grid(vtk_filename, vtk_points,
                                                                 vtk_cells)
  else
    vtk_nodedata = nothing
  end
  @timeit "build VTK grid (cell data)" vtk_celldata = vtk_grid(vtk_celldata_filename,
                                                                vtk_celldata_points,
                                                                vtk_celldata_cells)

  return vtk_nodedata, vtk_celldata
end


function calc_node_coordinates(mesh::TreeMesh, nodes, n_visnodes)
  coordinates, levels, _, _ = extract_mesh_information(mesh)

  # Extract number of spatial dimensions
  ndims_ = size(coordinates, 1)
  n_elements = length(levels)

  node_coordinates = Array{Float64, ndims_+2}(undef, ndims_, ntuple(_ -> n_visnodes, ndims_)..., n_elements)

  return calc_node_coordinates!(node_coordinates, nodes, mesh)
end


function calc_node_coordinates(mesh::StructuredMesh, nodes, n_visnodes)
  # Extract number of spatial dimensions
  ndims_ = ndims(mesh)
  n_elements = prod(size(mesh))

  f(args...; kwargs...) = Base.invokelatest(mesh.mapping, args...; kwargs...)

  node_coordinates = Array{Float64, ndims_+2}(undef, ndims_, ntuple(_ -> n_visnodes, ndims_)..., n_elements)

  return calc_node_coordinates!(node_coordinates, f, nodes, mesh)
end


function calc_node_coordinates(mesh::UnstructuredMesh2D, nodes, n_visnodes)
  # Extract number of spatial dimensions
  ndims_ = ndims(mesh)
  n_elements = length(mesh)

  # initialize the container for the node coordinates
  node_coordinates = Array{Float64, ndims_+2}(undef, ndims_, ntuple(_ -> n_visnodes, ndims_)..., n_elements)

  # work container for the corners of elements
  four_corners = zeros(eltype(mesh.corners), 4, 2)

  # loop through all elements and call the correct node coordinate constructor based on curved
  for element = 1:n_elements
    if mesh.element_is_curved[element]
      Trixi.calc_node_coordinates!(node_coordinates, element, nodes, view(mesh.surface_curves, :, element))
    else # straight sided element
      for i in 1:4, j in 1:2
        # pull the (x,y) values of these corners out of the global corners array
        four_corners[i, j] = mesh.corners[j, mesh.element_node_ids[i, element]]
      end
      Trixi.calc_node_coordinates!(node_coordinates, element, nodes, four_corners)
    end
  end

  return node_coordinates
end

# Version of calc_node_coordinates for a P4estMesh representing a manifold of dimension
# NDIMS embedded within an ambient space of dimension NDIMS_AMBIENT. This provides support 
# for standard 2D and 3D meshes (i.e. NDIMS = NDIMS_AMBIENT) as well as 2D surfaces in 3D 
# space (i.e. NDIMS = 2 and NDIMS_AMBIENT = 3).
function calc_node_coordinates(mesh::P4estMesh{NDIMS, NDIMS_AMBIENT}, nodes, 
                               n_visnodes) where {NDIMS, NDIMS_AMBIENT}

  node_coordinates = Array{Float64, NDIMS+2}(undef, NDIMS_AMBIENT,
                                              ntuple(_ -> n_visnodes, NDIMS)...,
                                              Trixi.ncells(mesh))

  return Trixi.calc_node_coordinates!(node_coordinates, mesh, nodes)
end


function calc_node_coordinates(mesh::T8codeMesh, nodes, n_visnodes)
  # Extract number of spatial dimensions
  ndims_ = ndims(mesh)

  node_coordinates = Array{Float64, ndims_+2}(undef, ndims_,
                                              ntuple(_ -> n_visnodes, ndims_)...,
                                              Trixi.ncells(mesh))

  return Trixi.calc_node_coordinates!(node_coordinates, mesh, nodes)
end


# Calculation of the node coordinates for `TreeMesh` in 2D
function calc_node_coordinates!(node_coordinates, nodes, mesh::TreeMesh{2})
  _, levels, _, length_level_0 = extract_mesh_information(mesh)

  # Extract number of spatial dimensions
  n_elements = length(levels)

  # Set the reference length
  # TODO: Is this ever different from 2?
  reference_length = 2.0
  # Compute the offset of the midpoint of the 1D reference interval
  # (its difference from zero)
  reference_offset = first(nodes) + reference_length / 2

  # Recompute the cell ids
  cell_ids = Trixi.local_leaf_cells(mesh.tree)

  # Calculate inverse Jacobian and node coordinates
  for element in 1:n_elements
    # Get cell id
    cell_id = cell_ids[element]

    # Get cell length
    dx = length_level_0 / 2^levels[element]

    # Get Jacobian value
    jacobian = dx / reference_length

    # Calculate node coordinates
    # Note that the `tree_coordinates` are the midpoints of the cells.
    # Hence, we need to add an offset for `nodes` with a midpoint
    # different from zero.
    for j in eachindex(nodes), i in eachindex(nodes)
      node_coordinates[1, i, j, element] = (
          mesh.tree.coordinates[1, cell_id] + jacobian * (nodes[i] - reference_offset))
      node_coordinates[2, i, j, element] = (
          mesh.tree.coordinates[2, cell_id] + jacobian * (nodes[j] - reference_offset))
    end
  end

  return node_coordinates
end


# Calculation of the node coordinates for `TreeMesh` in 3D
function calc_node_coordinates!(node_coordinates, nodes, mesh::TreeMesh{3})
  _, levels, _, length_level_0 = extract_mesh_information(mesh)

  # Extract number of spatial dimensions
  n_elements = length(levels)

  # Set the reference length
  # TODO: Is this ever different from 2?
  reference_length = 2.0
  # Compute the offset of the midpoint of the 1D reference interval
  # (its difference from zero)
  reference_offset = first(nodes) + reference_length / 2

  # Recompute the cell ids
  cell_ids = Trixi.local_leaf_cells(mesh.tree)

  # Calculate inverse Jacobian and node coordinates
  for element in 1:n_elements
    # Get cell id
    cell_id = cell_ids[element]

    # Get cell length
    dx = length_level_0 / 2^levels[element]

    # Get Jacobian value
    jacobian = dx / reference_length

    # Calculate node coordinates
    # Note that the `tree_coordinates` are the midpoints of the cells.
    # Hence, we need to add an offset for `nodes` with a midpoint
    # different from zero.
    for k in eachindex(nodes), j in eachindex(nodes), i in eachindex(nodes)
      node_coordinates[1, i, j, k, element] = (
          mesh.tree.coordinates[1, cell_id] + jacobian * (nodes[i] - reference_offset))
      node_coordinates[2, i, j, k, element] = (
          mesh.tree.coordinates[2, cell_id] + jacobian * (nodes[j] - reference_offset))
      node_coordinates[3, i, j, k, element] = (
          mesh.tree.coordinates[3, cell_id] + jacobian * (nodes[k] - reference_offset))
    end
  end

  return node_coordinates
end


# Calculation of the node coordinates for `StructuredMesh` in 2D
function calc_node_coordinates!(node_coordinates::AbstractArray{<:Any, 4}, f, nodes, mesh)
  linear_indices = LinearIndices(size(mesh))

  # Get cell length in reference mesh
  dx = 2 / size(mesh, 1)
  dy = 2 / size(mesh, 2)

  for cell_y in 1:size(mesh, 2), cell_x in 1:size(mesh, 1)
    element = linear_indices[cell_x, cell_y]

    # Calculate node coordinates of reference mesh
    cell_x_offset = -1 + (cell_x-1) * dx + dx/2
    cell_y_offset = -1 + (cell_y-1) * dy + dy/2

    for j in eachindex(nodes), i in eachindex(nodes)
      # node_coordinates are the mapped reference node_coordinates
      node_coordinates[:, i, j, element] .= f(cell_x_offset + dx/2 * nodes[i],
                                              cell_y_offset + dy/2 * nodes[j])
    end
  end

  return node_coordinates
end


# Calculation of the node coordinates for `StructuredMesh` in 3D
function calc_node_coordinates!(node_coordinates::AbstractArray{<:Any, 5}, f, nodes, mesh)
  linear_indices = LinearIndices(size(mesh))

  # Get cell length in reference mesh
  dx = 2 / size(mesh, 1)
  dy = 2 / size(mesh, 2)
  dz = 2 / size(mesh, 3)

  for cell_z in 1:size(mesh, 3), cell_y in 1:size(mesh, 2), cell_x in 1:size(mesh, 1)
    element = linear_indices[cell_x, cell_y, cell_z]

    # Calculate node coordinates of reference mesh
    cell_x_offset = -1 + (cell_x-1) * dx + dx/2
    cell_y_offset = -1 + (cell_y-1) * dy + dy/2
    cell_z_offset = -1 + (cell_z-1) * dz + dz/2

    for k in eachindex(nodes), j in eachindex(nodes), i in eachindex(nodes)
      # node_coordinates are the mapped reference node_coordinates
      node_coordinates[:, i, j, k, element] .= f(cell_x_offset + dx/2 * nodes[i],
                                                 cell_y_offset + dy/2 * nodes[j],
                                                 cell_z_offset + dz/2 * nodes[k])
    end
  end

  return node_coordinates
end


# Determine and return filenames for PVD fields
function pvd_filenames(filenames, pvd, output_directory)
  # Determine pvd filename
  if !isnothing(pvd)
    # Use filename if given on command line
    filename = pvd

    # Strip of directory/extension
    filename, _ = splitext(splitdir(filename)[2])
  else
    filename = get_pvd_filename(filenames)

    # If filename is empty, it means we were not able to determine an
    # appropriate file thus the user has to supply one
    if filename == ""
      error("could not auto-detect PVD filename (input file names have no common prefix): " *
            "please provide a PVD filename name with the keyword argument `pvd=...`")
    end
  end

  # Get full filenames
  pvd_filename = joinpath(output_directory, filename)
  pvd_celldata_filename = pvd_filename * "_celldata"

  return pvd_filename, pvd_celldata_filename
end


# Determine filename for PVD file based on common name
function get_pvd_filename(filenames)
  filenames = getindex.(splitdir.(filenames), 2)
  bases = getindex.(splitext.(filenames), 1)
  pvd_filename = longest_common_prefix(bases)
  return pvd_filename
end


# Determine longest common prefix
function longest_common_prefix(strings)
  # Return early if array is empty
  if isempty(strings)
    return ""
  end

  # Count length of common prefix, by ensuring that all strings are long enough
  # and then comparing the next character
  len = 0
  while all(length.(strings) .> len) && all(getindex.(strings, len+1) .== strings[1][len+1])
    len +=1
  end

  return strings[1][1:len]
end


# Convert coordinates and level information to a list of points and VTK cells (2D version)
function calc_vtk_points_cells(::Val{2}, coordinates::AbstractMatrix{Float64},
                               levels::AbstractVector{Int},
                               center_level_0::AbstractVector{Float64},
                               length_level_0::Float64,
                               n_visnodes::Int=1)
  ndim = 2
  # Create point locator
  pl = PointLocator{2}(center_level_0, length_level_0, 1e-12)

  # Create arrays for points and cells
  n_elements = length(levels)
  points = Vector{SVector{2, Float64}}()
  vtk_cells = Vector{MeshCell}(undef, n_elements * n_visnodes^ndim)
  point_ids = Vector{Int}(undef, 2^ndim)

  # Reshape cell array for easy-peasy access
  reshaped = reshape(vtk_cells, n_visnodes, n_visnodes, n_elements)

  # Create VTK cell for each Trixi element
  for element_id in 1:n_elements
    # Extract cell values
    cell_x = coordinates[1, element_id]
    cell_y = coordinates[2, element_id]
    cell_dx = length_level_0 / 2^levels[element_id]

    # Adapt to visualization nodes for easy-to-understand loops
    dx = cell_dx / n_visnodes
    x_lowerleft = cell_x - cell_dx/2 - dx/2
    y_lowerleft = cell_y - cell_dx/2 - dx/2

    # Create cell for each visualization node
    for j = 1:n_visnodes
      for i = 1:n_visnodes
        # Determine x and y
        x = x_lowerleft + i * dx
        y = y_lowerleft + j * dx

        # Get point id for each vertex
        point_ids[1] = insert!(pl, points, (x - dx/2, y - dx/2))
        point_ids[2] = insert!(pl, points, (x + dx/2, y - dx/2))
        point_ids[3] = insert!(pl, points, (x - dx/2, y + dx/2))
        point_ids[4] = insert!(pl, points, (x + dx/2, y + dx/2))

        # Add cell
        reshaped[i, j, element_id] = MeshCell(VTKCellTypes.VTK_PIXEL, copy(point_ids))
      end
    end
  end

  # Convert array-of-points to two-dimensional array
  vtk_points = Matrix{Float64}(undef, ndim, length(points))
  for point_id in 1:length(points)
    vtk_points[1, point_id] = points[point_id][1]
    vtk_points[2, point_id] = points[point_id][2]
  end

  return vtk_points, vtk_cells
end


# Convert coordinates and level information to a list of points and VTK cells (3D version)
function calc_vtk_points_cells(::Val{3}, coordinates::AbstractMatrix{Float64},
                               levels::AbstractVector{Int},
                               center_level_0::AbstractVector{Float64},
                               length_level_0::Float64,
                               n_visnodes::Int=1)
  ndim = 3
  # Create point locator
  pl = PointLocator{3}(center_level_0, length_level_0, 1e-12)

  # Create arrays for points and cells
  n_elements = length(levels)
  points = Vector{SVector{3, Float64}}()
  vtk_cells = Vector{MeshCell}(undef, n_elements * n_visnodes^ndim)
  point_ids = Vector{Int}(undef, 2^ndim)

  # Reshape cell array for easy-peasy access
  reshaped = reshape(vtk_cells, n_visnodes, n_visnodes, n_visnodes, n_elements)

  # Create VTK cell for each Trixi element
  for element_id in 1:n_elements
    # Extract cell values
    cell_x = coordinates[1, element_id]
    cell_y = coordinates[2, element_id]
    cell_z = coordinates[3, element_id]
    cell_dx = length_level_0 / 2^levels[element_id]

    # Adapt to visualization nodes for easy-to-understand loops
    dx = cell_dx / n_visnodes
    x_lowerleft = cell_x - cell_dx/2 - dx/2
    y_lowerleft = cell_y - cell_dx/2 - dx/2
    z_lowerleft = cell_z - cell_dx/2 - dx/2

    # Create cell for each visualization node
    for k = 1:n_visnodes, j = 1:n_visnodes, i = 1:n_visnodes
      # Determine x and y
      x = x_lowerleft + i * dx
      y = y_lowerleft + j * dx
      z = z_lowerleft + k * dx

      # Get point id for each vertex
      point_ids[1] = insert!(pl, points, (x - dx/2, y - dx/2, z - dx/2))
      point_ids[2] = insert!(pl, points, (x + dx/2, y - dx/2, z - dx/2))
      point_ids[3] = insert!(pl, points, (x - dx/2, y + dx/2, z - dx/2))
      point_ids[4] = insert!(pl, points, (x + dx/2, y + dx/2, z - dx/2))
      point_ids[5] = insert!(pl, points, (x - dx/2, y - dx/2, z + dx/2))
      point_ids[6] = insert!(pl, points, (x + dx/2, y - dx/2, z + dx/2))
      point_ids[7] = insert!(pl, points, (x - dx/2, y + dx/2, z + dx/2))
      point_ids[8] = insert!(pl, points, (x + dx/2, y + dx/2, z + dx/2))

      # Add cell
      reshaped[i, j, k, element_id] = MeshCell(VTKCellTypes.VTK_VOXEL, copy(point_ids))
    end
  end

  # Convert array-of-points to two-dimensional array
  vtk_points = Matrix{Float64}(undef, ndim, length(points))
  for point_id in 1:length(points)
    vtk_points[1, point_id] = points[point_id][1]
    vtk_points[2, point_id] = points[point_id][2]
    vtk_points[3, point_id] = points[point_id][3]
  end

  return vtk_points, vtk_cells
end


# Convert coordinates and level information to a list of points and VTK cells for `StructuredMesh` (2D version)
function calc_vtk_points_cells(node_coordinates::AbstractArray{<:Any,4})
  n_elements = size(node_coordinates, 4)
  size_ = size(node_coordinates)
  n_points = prod(size_[2:end])
  # Linear indices to access points by node indices and element id
  linear_indices = LinearIndices(size_[2:end])

  # Use Lagrange nodes as VTK points. Note that we call size(node_coordinates, 1) in order  
  # to provide support for standard two-dimensional meshes as well as meshes representing 
  # 2D surfaces in 3D space, which are implemented using P4estMesh{2, 3}.
  vtk_points = reshape(node_coordinates, (size(node_coordinates, 1), n_points))
  vtk_cells = Vector{MeshCell}(undef, n_elements)

  # Create cell for each element
  for element in 1:n_elements
    vertices = [linear_indices[1, 1, element],
                linear_indices[end, 1, element],
                linear_indices[end, end, element],
                linear_indices[1, end, element]]

    edges = vcat(linear_indices[2:end-1, 1, element],
                 linear_indices[end, 2:end-1, element],
                 linear_indices[2:end-1, end, element],
                 linear_indices[1, 2:end-1, element])

    faces = vec(linear_indices[2:end-1, 2:end-1, element])

    point_ids = vcat(vertices, edges, faces)
    vtk_cells[element] = MeshCell(VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL, point_ids)
  end

  return vtk_points, vtk_cells
end


# Convert coordinates and level information to a list of points and VTK cells for `StructuredMesh` (3D version)
function calc_vtk_points_cells(node_coordinates::AbstractArray{<:Any,5})
  n_elements = size(node_coordinates, 5)
  size_ = size(node_coordinates)
  n_points = prod(size_[2:end])
  # Linear indices to access points by node indices and element id
  linear_indices = LinearIndices(size_[2:end])

  # Use lagrange nodes as VTK points
  vtk_points = reshape(node_coordinates, (3, n_points))
  vtk_cells = Vector{MeshCell}(undef, n_elements)

  # Create cell for each element
  for element in 1:n_elements
    vertices = [linear_indices[1, 1, 1, element],
                linear_indices[end, 1, 1, element],
                linear_indices[end, end, 1, element],
                linear_indices[1, end, 1, element],
                linear_indices[1, 1, end, element],
                linear_indices[end, 1, end, element],
                linear_indices[end, end, end, element],
                linear_indices[1, end, end, element]]

    # This order doesn't make any sense. This is completely different
    # from what is shown in
    # https://blog.kitware.com/wp-content/uploads/2018/09/Source_Issue_43.pdf
    # but this is the way it works.
    edges = vcat(linear_indices[2:end-1, 1, 1, element],
                 linear_indices[end, 2:end-1, 1, element],
                 linear_indices[2:end-1, end, 1, element],
                 linear_indices[1, 2:end-1, 1, element],
                 linear_indices[2:end-1, 1, end, element],
                 linear_indices[end, 2:end-1, end, element],
                 linear_indices[2:end-1, end, end, element],
                 linear_indices[1, 2:end-1, end, element],
                 linear_indices[1, 1, 2:end-1, element],
                 linear_indices[end, 1, 2:end-1, element],
                 linear_indices[1, end, 2:end-1, element],
                 linear_indices[end, end, 2:end-1, element])

    # See above
    faces = vcat(vec(linear_indices[1, 2:end-1, 2:end-1, element]),
                 vec(linear_indices[end, 2:end-1, 2:end-1, element]),
                 vec(linear_indices[2:end-1, 1, 2:end-1, element]),
                 vec(linear_indices[2:end-1, end, 2:end-1, element]),
                 vec(linear_indices[2:end-1, 2:end-1, 1, element]),
                 vec(linear_indices[2:end-1, 2:end-1, end, element]))

    volume = vec(linear_indices[2:end-1, 2:end-1, 2:end-1, element])

    point_ids = vcat(vertices, edges, faces, volume)
    vtk_cells[element] = MeshCell(VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON, point_ids)
  end

  return vtk_points, vtk_cells
end

# Check if file is a data file
function is_solution_restart_file(filename::String)
  # Open file for reading
  h5open(filename, "r") do file
    # If attribute "mesh_file" exists, this must be a data file
    return haskey(attributes(file), "mesh_file")
  end
end


# Use data file to extract mesh filename from attributes
function extract_mesh_filename(filename::String)
  # Open file for reading
  h5open(filename, "r") do file
    # Extract filename relative to data file
    mesh_file = read(attributes(file)["mesh_file"])

    return joinpath(dirname(filename), mesh_file)
  end
end


function extract_mesh_information(mesh::TreeMesh)
  center_level_0 = mesh.tree.center_level_0
  length_level_0 = mesh.tree.length_level_0

  coordinates = mesh.tree.coordinates[:, Trixi.leaf_cells(mesh.tree)]
  levels = mesh.tree.levels[Trixi.leaf_cells(mesh.tree)]

  return coordinates, levels, center_level_0, length_level_0
end


# Read in data file and return all relevant information
function read_datafile(filename::String, ::Union{TreeMesh, StructuredMesh,
                                                 UnstructuredMesh2D, P4estMesh, T8codeMesh})
  # Open file for reading
  h5open(filename, "r") do file
    # Extract basic information
    ndims_ = read(attributes(file)["ndims"])
    polydeg = read(attributes(file)["polydeg"])
    n_elements = read(attributes(file)["n_elements"])
    n_variables = read(attributes(file)["n_vars"])
    time = read(attributes(file)["time"])

    # Extract labels for legend
    labels = Array{String}(undef, 1, n_variables)
    for v = 1:n_variables
      labels[1, v] = read(attributes(file["variables_$v"])["name"])
    end

    # Extract data arrays
    n_nodes = polydeg + 1

    if ndims_ == 2
      data = Array{Float64}(undef, n_nodes, n_nodes, n_elements, n_variables)
      for v = 1:n_variables
        vardata = read(file["variables_$v"])
        @views data[:, :, :, v][:] .= vardata
      end
    elseif ndims_ == 3
      data = Array{Float64}(undef, n_nodes, n_nodes, n_nodes, n_elements, n_variables)
      for v = 1:n_variables
        vardata = read(file["variables_$v"])
        @views data[:, :, :, :, v][:] .= vardata
      end
    else
      error("Unsupported number of spatial dimensions: ", ndims_)
    end

    # Extract element variable arrays
    element_variables = Dict{String, Union{Vector{Float64}, Vector{Int}}}()
    index = 1
    while haskey(file, "element_variables_$index")
      varname = read(attributes(file["element_variables_$index"])["name"])
      element_variables[varname] = read(file["element_variables_$index"])
      index +=1
    end

    # Extract node variable arrays
    node_variables = Dict{String, Union{Array{Float64}, Array{Int}}}()
    index = 1
    while haskey(file, "node_variables_$index")
      varname = read(attributes(file["node_variables_$index"])["name"])
      nodedata = read(file["node_variables_$index"])
      if ndims_ == 2
        node_variables[varname] = Array{Float64}(undef, n_nodes, n_nodes, n_elements)
        @views node_variables[varname][:, :, :] .= nodedata
      else
        error("Unsupported number of spatial dimensions: ", ndims_)
      end
      index +=1
    end

    return labels, data, n_elements, n_nodes, element_variables, node_variables, time
  end
end

# Read in data file and return all relevant information for a DGMulti solution
function read_datafile(filename::String, ::DGMultiMesh)
  # Open file for reading
  h5open(filename, "r") do file
    # Extract basic information
    ndims_ = read(attributes(file)["ndims"])

    etype_str = read(attributes(file)["element_type"])
    etype = Trixi.get_element_type_from_string(etype_str)()

    if etype isa Trixi.Wedge && haskey(attributes(file), "polydeg_tri")
        polydeg = tuple(read(attributes(file)["polydeg_tri"]),
                        read(attributes(file)["polydeg_line"]))
    else
        polydeg = read(attributes(file)["polydeg"])
    end

    n_elements = read(attributes(file)["n_elements"])
    dof_per_elem = read(attributes(file)["dof_per_elem"])
    n_variables = read(attributes(file)["n_vars"])
    time = read(attributes(file)["time"])

    # Extract labels for legend
    labels = Array{String}(undef, 1, n_variables)
    for v = 1:n_variables
      labels[1, v] = read(attributes(file["variables_$v"])["name"])
    end

    # Extract data arrays
    n_nodes = maximum(polydeg) + 1

    data = Array{Float64}(undef, dof_per_elem, n_elements, n_variables)
    for v = 1:n_variables
      vardata = read(file["variables_$v"])
      @views data[:, :, v] .= vardata
    end

    # Extract element variable arrays
    element_variables = Dict{String, Union{Vector{Float64}, Vector{Int}}}()
    index = 1
    while haskey(file, "element_variables_$index")
      varname = read(attributes(file["element_variables_$index"])["name"])
      element_variables[varname] = read(file["element_variables_$index"])
      index +=1
    end

    # Extract node variable arrays
    node_variables = Dict{String, Union{Array{Float64}, Array{Int}}}()
    index = 1
    while haskey(file, "node_variables_$index")
      varname = read(attributes(file["node_variables_$index"])["name"])
      nodedata = read(file["node_variables_$index"])
      if ndims_ == 2
        node_variables[varname] = Array{Float64}(undef, n_nodes, n_nodes, n_elements)
        @views node_variables[varname][:, :, :] .= nodedata
      else
        error("Unsupported number of spatial dimensions: ", ndims_)
      end
      index +=1
    end

    return labels, data, n_elements, n_nodes, element_variables, node_variables, time, etype, polydeg
  end
end

# Reconstruct a `RefElemData` from information in `mesh_file`
function load_basis(mesh_file, ::DGMultiMesh)
  etype_str = h5open(mesh_file, "r") do file
      return read(attributes(file)["element_type"])
  end

  etype = Trixi.get_element_type_from_string(etype_str)()
  polydeg = h5open(mesh_file, "r") do file
      if etype isa Trixi.Wedge && haskey(attributes(file), "polydeg_tri")
          return tuple(read(attributes(file)["polydeg_tri"]),
                      read(attributes(file)["polydeg_line"]))
      else
          return read(attributes(file)["polydeg"])
      end
  end

  if etype isa Trixi.Wedge && polydeg isa NTuple{2}
    factor_a = Trixi.RefElemData(Trixi.StartUpDG.Tri(), Trixi.Polynomial(), polydeg[1])
    factor_b = Trixi.RefElemData(Trixi.StartUpDG.Line(), Trixi.Polynomial(), polydeg[2])

    tensor = Trixi.TensorProductWedge(factor_a, factor_b)
    rd = Trixi.RefElemData(etype, tensor)
  else
    rd = Trixi.RefElemData(etype, Trixi.Polynomial(), polydeg)
  end

  return rd
end

# For other mesh types, return nothing
function load_basis(mesh_file, ::Union{TreeMesh, StructuredMesh,
                                       UnstructuredMesh2D, P4estMesh, 
                                       T8codeMesh})
  return nothing
end
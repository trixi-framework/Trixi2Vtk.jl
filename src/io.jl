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
function read_datafile(filename::String)
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

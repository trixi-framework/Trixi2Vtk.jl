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


# Read in mesh file and return relevant data
function read_meshfile(filename::String)
  # Open file for reading
  h5open(filename, "r") do file
    # Extract basic information
    if haskey(attributes(file), "ndims")
      ndims_ = read(attributes(file)["ndims"])
    else
      ndims_ = read(attributes(file)["ndim"]) # FIXME once Trixi's 3D branch is merged & released
    end
    n_cells = read(attributes(file)["n_cells"])
    n_leaf_cells = read(attributes(file)["n_leaf_cells"])
    center_level_0 = read(attributes(file)["center_level_0"])
    length_level_0 = read(attributes(file)["length_level_0"])

    # Extract coordinates, levels, child cells
    coordinates = Array{Float64}(undef, ndims_, n_cells)
    coordinates .= read(file["coordinates"])
    levels = Array{Int}(undef, n_cells)
    levels .= read(file["levels"])
    child_ids = Array{Int}(undef, 2^ndims_, n_cells)
    child_ids .= read(file["child_ids"])

    # Extract leaf cells (= cells to be plotted) and contract all other arrays accordingly
    leaf_cells = similar(levels)
    n_cells = 0
    for cell_id in 1:length(levels)
      if sum(child_ids[:, cell_id]) > 0
        continue
      end

      n_cells += 1
      leaf_cells[n_cells] = cell_id
    end
    leaf_cells = leaf_cells[1:n_cells]

    coordinates = coordinates[:, leaf_cells]
    levels = levels[leaf_cells]

    return center_level_0, length_level_0, leaf_cells, coordinates, levels
  end
end


# Read in structured mesh file and return relevant data
function read_meshfile_structured(filename::String; RealT=Float64)
  # Open file for reading
  h5open(filename, "r") do file
    size_ = Tuple(read(attributes(file)["size"]))
    mapping_as_string = read(attributes(file)["mapping"])

    # A temporary workaround to evaluate the code that defines the domain mapping in a local scope.
    # This prevents errors when multiple restart elixirs are executed in one session, where one
    # defines `mapping` as a variable, while the other defines it as a function.
    #
    # This should be replaced with something more robust and secure,
    # see https://github.com/trixi-framework/Trixi.jl/issues/541).
    expr = Meta.parse(mapping_as_string)
    if expr.head == :toplevel
      expr.head = :block
    end

    mapping = @eval function(xi, eta)
      $expr
      mapping(xi, eta)
    end

    return Trixi.CurvedMesh(size_, mapping; RealT=RealT, unsaved_changes=false,
                            mapping_as_string=mapping_as_string)
  end
end


# Read in data file and return all relevant information
function read_datafile(filename::String)
  # Open file for reading
  h5open(filename, "r") do file
    # Extract basic information
    if haskey(attributes(file), "ndims")
      ndims_ = read(attributes(file)["ndims"])
    else
      ndims_ = read(attributes(file)["ndim"])
    end
    if haskey(attributes(file), "polydeg")
      polydeg = read(attributes(file)["polydeg"])
    else
      polydeg = read(attributes(file)["N"])
    end
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

    return labels, data, n_elements, n_nodes, element_variables, time
  end
end


# Read in data file and return all relevant information
function read_datafile_structured(filename::String)
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

    return labels, data, n_elements, n_nodes, element_variables, time
  end
end


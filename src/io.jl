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


# TODO Trixi needs to be adapted to work with mesh files instead of restart files.
# Then, these functions can be deleted and the ones from Trixi can be used instead.
function load_mesh_serial(filename::AbstractString; RealT, n_cells_max=0)
  ndims, mesh_type = h5open(filename, "r") do file
    return read(attributes(file)["ndims"]),
           read(attributes(file)["mesh_type"])
  end

  if mesh_type == "TreeMesh"
    n_cells = h5open(filename, "r") do file
      return read(attributes(file)["n_cells"])
    end
    n_cells_max = max(n_cells_max, n_cells)
    mesh = Trixi.TreeMesh(Trixi.SerialTree{ndims}, n_cells_max)
    load_mesh!(mesh, filename)
  elseif mesh_type == "CurvedMesh"
    size_, mapping_as_string = h5open(filename, "r") do file
      return read(attributes(file)["size"]),
             read(attributes(file)["mapping"])
    end

    size = Tuple(size_)

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

    if ndims == 1
      mapping = @eval function(xi)
        $expr
        mapping(xi)
      end
    elseif ndims == 2
      mapping = @eval function(xi, eta)
        $expr
        mapping(xi, eta)
      end
    else # ndims == 3
      mapping = @eval function(xi, eta, zeta)
        $expr
        mapping(xi, eta, zeta)
      end
    end

    mesh = Trixi.CurvedMesh(size, mapping; RealT=RealT, unsaved_changes=false, mapping_as_string=mapping_as_string)
  else
    error("Unknown mesh type!")
  end

  return mesh
end


function load_mesh!(mesh::Trixi.TreeMesh, filename::AbstractString)
  # Determine mesh filename
  mesh.current_filename = filename
  mesh.unsaved_changes = false

  # Read mesh file
  h5open(filename, "r") do file
    # Set domain information
    mesh.tree.center_level_0 = read(attributes(file)["center_level_0"])
    mesh.tree.length_level_0 = read(attributes(file)["length_level_0"])
    mesh.tree.periodicity    = Tuple(read(attributes(file)["periodicity"]))

    # Set length
    n_cells = read(attributes(file)["n_cells"])
    resize!(mesh.tree, n_cells)

    # Read in data
    mesh.tree.parent_ids[1:n_cells] = read(file["parent_ids"])
    mesh.tree.child_ids[:, 1:n_cells] = read(file["child_ids"])
    mesh.tree.neighbor_ids[:, 1:n_cells] = read(file["neighbor_ids"])
    mesh.tree.levels[1:n_cells] = read(file["levels"])
    mesh.tree.coordinates[:, 1:n_cells] = read(file["coordinates"])
  end

  return mesh
end


function extract_mesh_information(mesh::Trixi.TreeMesh)
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

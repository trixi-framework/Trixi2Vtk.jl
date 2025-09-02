
# Interpolate data from input format to desired output format (vtu version)
function interpolate_data(::Val{:vtu}, input_data, mesh::TreeMesh, n_visnodes, verbose)
  # Calculate equidistant output nodes (visualization nodes)
  dx = 2 / n_visnodes
  nodes_out = collect(range(-1 + dx/2, 1 - dx/2, length=n_visnodes))

  return raw2interpolated(input_data, nodes_out)
end


# Interpolate data from input format to desired output format (StructuredMesh or UnstructuredMesh2D version)
function interpolate_data(::Val{:vtu}, input_data,
                          mesh::Union{StructuredMesh, UnstructuredMesh2D, P4estMesh, T8codeMesh},
                          n_visnodes, verbose)
  # Calculate equidistant output nodes
  nodes_out = collect(range(-1, 1, length=n_visnodes))

  return raw2interpolated(input_data, nodes_out)
end


# Interpolate data from input format to desired output format (vti version)
function interpolate_data(::Val{:vti}, input_data, mesh::TreeMesh, n_visnodes, verbose)
  coordinates, levels, center_level_0, length_level_0 = extract_mesh_information(mesh)

  # Normalize element coordinates: move center to origin and domain size to [-1, 1]²
  normalized_coordinates = similar(coordinates)
  for element_id in axes(coordinates, 2)
    @views normalized_coordinates[:, element_id] .= (
        (coordinates[:, element_id] .- center_level_0) ./ (length_level_0 / 2 ))
  end

  # Determine level-wise resolution
  max_level = maximum(levels)
  resolution = n_visnodes * 2^max_level

  # nvisnodes_per_level is an array (accessed by "level + 1" to accommodate
  # level-0-cell) that contains the number of visualization nodes for any
  # refinement level to visualize on an equidistant grid
  nvisnodes_per_level = [2^(max_level - level)*n_visnodes for level in 0:max_level]

  # Interpolate unstructured DG data to structured data
  structured_data = unstructured2structured(input_data, normalized_coordinates, levels,
                                            resolution, nvisnodes_per_level)

  return structured_data
end


# Interpolate unstructured DG data to structured data (cell-centered)
function unstructured2structured(unstructured_data::AbstractArray{Float64},
                                 normalized_coordinates::AbstractArray{Float64},
                                 levels::AbstractArray{Int}, resolution::Int,
                                 nvisnodes_per_level::AbstractArray{Int})
  # Extract number of spatial dimensions
  ndims_ = size(normalized_coordinates, 1)

  # Extract data shape information
  n_nodes_in, _, n_elements, n_variables = size(unstructured_data)

  # Get node coordinates for DG locations on reference element
  nodes_in, _ = gauss_lobatto_nodes_weights(n_nodes_in)

  # Calculate interpolation vandermonde matrices for each level
  max_level = length(nvisnodes_per_level) - 1
  vandermonde_per_level = []
  for l in 0:max_level
    n_nodes_out = nvisnodes_per_level[l + 1]
    dx = 2 / n_nodes_out
    nodes_out = collect(range(-1 + dx/2, 1 - dx/2, length=n_nodes_out))
    push!(vandermonde_per_level, polynomial_interpolation_matrix(nodes_in, nodes_out))
  end

  # For each element, calculate index position at which to insert data in global data structure
  lower_left_index = element2index(normalized_coordinates, levels, resolution, nvisnodes_per_level)

  # Create output data structure
  structured = Array{Float64}(undef, resolution, resolution, n_variables)
  # For each variable, interpolate element data and store to global data structure
  for v in 1:n_variables
    # Reshape data array for use in interpolate_nodes function
    reshaped_data = reshape(unstructured_data[:, :, :, v], 1, n_nodes_in, n_nodes_in, n_elements)

    for element_id in 1:n_elements
      # Extract level for convenience
      level = levels[element_id]

      # Determine target indices
      n_nodes_out = nvisnodes_per_level[level + 1]
      first = lower_left_index[:, element_id]
      last = first .+ (n_nodes_out - 1)

      # Interpolate data
      vandermonde = vandermonde_per_level[level + 1]
      structured[first[1]:last[1], first[2]:last[2], v] .= (
          reshape(interpolate_nodes(reshaped_data[:, :, :, element_id], vandermonde, 1),
                  n_nodes_out, n_nodes_out))
    end
  end

  # Return as one 1D array for each variable
  return reshape(structured, resolution^ndims_, n_variables)
end


# For a given normalized element coordinate, return the index of its lower left
# contribution to the global data structure
function element2index(normalized_coordinates::AbstractArray{Float64}, levels::AbstractArray{Int},
                       resolution::Int, nvisnodes_per_level::AbstractArray{Int})
  # Extract number of spatial dimensions
  ndims_ = size(normalized_coordinates, 1)

  n_elements = length(levels)

  # First, determine lower left coordinate for all cells
  dx = 2 / resolution
  lower_left_coordinate = Array{Float64}(undef, ndims_, n_elements)
  for element_id in 1:n_elements
    nvisnodes = nvisnodes_per_level[levels[element_id] + 1]
    lower_left_coordinate[1, element_id] = (
        normalized_coordinates[1, element_id] - (nvisnodes - 1)/2 * dx)
    lower_left_coordinate[2, element_id] = (
        normalized_coordinates[2, element_id] - (nvisnodes - 1)/2 * dx)
  end

  # Then, convert coordinate to global index
  indices = coordinate2index(lower_left_coordinate, resolution)

  return indices
end


# Find 2D array index for a 2-tuple of normalized, cell-centered coordinates (i.e., in [-1,1])
function coordinate2index(coordinate, resolution::Integer)
  # Calculate 1D normalized coordinates
  dx = 2 / resolution
  mesh_coordinates = collect(range(-1 + dx/2, 1 - dx/2, length=resolution))

  # Find index
  id_x = searchsortedfirst.(Ref(mesh_coordinates), coordinate[1, :], lt=(x,y)->x .< y .- dx/2)
  id_y = searchsortedfirst.(Ref(mesh_coordinates), coordinate[2, :], lt=(x,y)->x .< y .- dx/2)
  return transpose(hcat(id_x, id_y))
end


# Interpolate to specified output nodes
function raw2interpolated(data_gl::AbstractArray{Float64}, nodes_out)
  # Extract number of spatial dimensions
  ndims_ = ndims(data_gl) - 2

  # Extract data shape information
  n_nodes_in = size(data_gl, 1)
  n_nodes_out = length(nodes_out)
  n_elements = size(data_gl, ndims_ + 1)
  n_variables = size(data_gl, ndims_ + 2)

  # Get node coordinates for DG locations on reference element
  nodes_in, _ = gauss_lobatto_nodes_weights(n_nodes_in)

  # Calculate Vandermonde matrix
  vandermonde = polynomial_interpolation_matrix(nodes_in, nodes_out)

  @assert ndims_ == 2 || ndims_ == 3 "Only 2D and 3D data supported"
  # Create output data structure
  data_vis = Array{Float64}(undef, ntuple(_ -> n_nodes_out, ndims_)..., n_elements, n_variables)

  # For each variable, interpolate element data and store to global data structure
  for v in 1:n_variables
    # Reshape data array for use in interpolate_nodes function
    @views reshaped_data = reshape(data_gl[.., v], 1,
                                   ntuple(_ -> n_nodes_in, ndims_)..., n_elements)

    # Interpolate data to visualization nodes
    for element_id in 1:n_elements
      @views data_vis[.., element_id, v] .= reshape(
          interpolate_nodes(reshaped_data[.., element_id], vandermonde, 1),
          ntuple(_ -> n_nodes_out, ndims_)...)
    end
  end

  # Return as one 1D array for each variable
  return reshape(data_vis, n_nodes_out^ndims_ * n_elements, n_variables)
end

# Interpolate data using the given Vandermonde matrix and return interpolated values (2D version).
function interpolate_nodes(data_in::AbstractArray{T, 3},
                           vandermonde, n_vars) where T
  n_nodes_out = size(vandermonde, 1)
  data_out = zeros(eltype(data_in), n_vars, n_nodes_out, n_nodes_out)
  interpolate_nodes!(data_out, data_in, vandermonde, n_vars)
end


function interpolate_nodes!(data_out::AbstractArray{T, 3}, data_in::AbstractArray{T, 3},
                            vandermonde, n_vars) where T
  n_nodes_out = size(vandermonde, 1)
  n_nodes_in  = size(vandermonde, 2)

  for j in 1:n_nodes_out
    for i in 1:n_nodes_out
      for v in 1:n_vars
        acc = zero(eltype(data_out))
        for jj in 1:n_nodes_in
          for ii in 1:n_nodes_in
            acc += vandermonde[i, ii] * data_in[v, ii, jj] * vandermonde[j, jj]
          end
        end
        data_out[v, i, j] = acc
      end
    end
  end

  return data_out
end


# Interpolate data using the given Vandermonde matrix and return interpolated values (3D version).
function interpolate_nodes(data_in::AbstractArray{T, 4},
                           vandermonde, n_vars) where T
  n_nodes_out = size(vandermonde, 1)
  data_out = zeros(eltype(data_in), n_vars, n_nodes_out, n_nodes_out, n_nodes_out)
  interpolate_nodes!(data_out, data_in, vandermonde, n_vars)
end

function interpolate_nodes!(data_out::AbstractArray{T, 4}, data_in::AbstractArray{T, 4},
                            vandermonde, n_vars) where T
  n_nodes_out = size(vandermonde, 1)
  n_nodes_in  = size(vandermonde, 2)

  for k in 1:n_nodes_out, j in 1:n_nodes_out, i in 1:n_nodes_out
    for v in 1:n_vars
      acc = zero(eltype(data_out))
      for kk in 1:n_nodes_in, jj in 1:n_nodes_in, ii in 1:n_nodes_in
        acc += vandermonde[i, ii] * vandermonde[j, jj] * vandermonde[k, kk] * data_in[v, ii, jj, kk]
      end
      data_out[v, i, j, k] = acc
    end
  end

  return data_out
end

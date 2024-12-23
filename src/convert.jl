"""
    trixi2vtk(filename::AbstractString...;
              format=:vtu, verbose=false, hide_progress=false, pvd=nothing,
              output_directory=".", nvisnodes=nothing, save_celldata=true,
              reinterpolate=true, data_is_uniform=false)

Convert Trixi-generated output files to VTK files (VTU or VTI).

# Arguments
- `filename`: One or more Trixi solution/restart/mesh files to convert to a VTK file.
              Filenames support file globbing, e.g., "solution*" to match all files starting
              with `solution`.
- `format`: Output format for solution/restart files. Can be 'vtu' or 'vti'.
- `verbose`: Set to `true` to enable verbose output.
- `hide_progress`: Hide progress bar (will be hidden automatically if `verbose` is `true`).
- `pvd`: Use this filename to store PVD file (instead of auto-detecting name). Note that
         only the name will be used (directory and file extension are ignored).
- `output_directory`: Output directory where generated files are stored.
- `nvisnodes`: Number of visualization nodes per element.
               (default: number of DG nodes for StructuredMesh or UnstructuredMesh2D,
                         twice the number of DG nodes for TreeMesh).
               A value of `0` (zero) uses the number of nodes in the DG elements.
- `save_celldata`: Boolean value to determine if cell-based data should be saved.
                   (default: `true`)
- `reinterpolate`: Boolean value to determine if data should be reinterpolated
                   onto uniform points. When `false` the raw data at the compute nodes
                   is copied into the appropriate format.
                   (default: `true`)
- `data_is_uniform`: Boolean to indicate if the data to be converted is from a finite difference
                     method on a uniform grid of points.
                     (default: `false`)

# Examples
```julia
julia> trixi2vtk("out/solution_000*.h5")
[...]
```
"""
function trixi2vtk(filename::AbstractString...;
                   format=:vtu, verbose=false, hide_progress=false, pvd=nothing,
                   output_directory=".", nvisnodes=nothing, save_celldata=true,
                   reinterpolate=true, data_is_uniform=false)
  # Reset timer
  reset_timer!()

  # Convert filenames to a single list of strings
  if isempty(filename)
    error("no input file was provided")
  end
  filenames = expand_filename_patterns(filename)
  if isempty(filenames)
    error("no such file(s): ", join(filename, ", "))
  end

  # Ensure valid format
  if !(format in (:vtu, :vti))
    error("unsupported output format '$format' (must be 'vtu' or 'vti')")
  end

  # If verbose mode is enabled, always hide progress bar
  if verbose
    hide_progress = true
  end

  # Variable to avoid writing PVD files if only a single file is converted
  is_single_file = length(filenames) == 1

  # Get pvd filenames and open files
  if !is_single_file
    pvd_filename, pvd_celldata_filename = pvd_filenames(filenames, pvd, output_directory)
    verbose && println("Opening PVD files '$(pvd_filename).pvd' + '$(pvd_celldata_filename).pvd'...")
    @timeit "open PVD file" begin
      pvd = paraview_collection(pvd_filename)
      pvd_celldata = paraview_collection(pvd_celldata_filename)
    end
  end

  # Variable to avoid writing PVD file if only mesh files were converted
  has_data = false

  # Show progress bar if not disabled
  if !hide_progress
    progress = Progress(length(filenames);
                        dt = 0.5,
                        desc = "Converting .h5 to .$(format)...",
                        barlen = 40)
  end

  # Show warning when reinterpolating node-level data of subcell limiting
  # Auxiliary variable to show warning only once
  has_warned_about_interpolation = false

  # Iterate over input files
  for (index, filename) in enumerate(filenames)
    verbose && println("Processing file $filename ($(index)/$(length(filenames)))...")

    # Check if data file exists
    if !isfile(filename)
      error("data file '$filename' does not exist")
    end

    # Check if it is a data file at all
    is_datafile = is_solution_restart_file(filename)

    # If file is solution/restart file, extract mesh file name
    if is_datafile
      # Get mesh file name
      meshfile = extract_mesh_filename(filename)

      # Check if mesh file exists
      if !isfile(meshfile)
        error("mesh file '$meshfile' does not exist")
      end
    else
      meshfile = filename
    end

    # Read mesh
    verbose && println("| Reading mesh file...")
    @timeit "read mesh" mesh = Trixi.load_mesh_serial(meshfile; n_cells_max=0, RealT=Float64)

    # Check compatibility of the mesh type and the output format
    if format === :vti && !(mesh isa Trixi.TreeMesh{2})
      throw(ArgumentError("VTI format only available for 2D TreeMesh"))
    end

    # Read data only if it is a data file
    if is_datafile
      verbose && println("| Reading data file...")
      @timeit "read data" (labels, data, n_elements, n_nodes,
                           element_variables, node_variables, time) = read_datafile(filename)

      assert_cells_elements(n_elements, mesh, filename, meshfile)

      # Determine resolution for data interpolation
      n_visnodes = get_default_nvisnodes_solution(nvisnodes, n_nodes, mesh)

      # If a user requests that no reinterpolation is done automatically set
      # `n_visnodes` to be the same as the number of nodes in the raw data.
      if !reinterpolate
        n_visnodes = n_nodes
      end

      # Check if the raw data is uniform (finite difference) or not (dg)
      # and create the corresponding node set for reinterpolation / copying.
      if (reinterpolate && !data_is_uniform) || (!reinterpolate && data_is_uniform)
        # (1) Default settings; presumably the most common
        # (2) Finite difference data
        node_set = collect(range(-1, 1, length=n_visnodes))
      elseif !reinterpolate && !data_is_uniform
        # raw data is on a set of LGL nodes
        node_set, _ = gauss_lobatto_nodes_weights(n_visnodes)
      else # reinterpolate & data_is_uniform
        throw(ArgumentError("Uniform data should not be reinterpolated! Set `reinterpolate=false` and try again."))
      end
    else
      # If file is a mesh file, do not interpolate data as detailed
      n_visnodes = get_default_nvisnodes_mesh(nvisnodes, mesh)
      # Create an "empty" node set that is unused in the mesh conversion
      node_set = Array{Float64}(undef, n_visnodes)
    end

    # Create output directory if it does not exist
    mkpath(output_directory)

    # Build VTK grids
    vtk_nodedata, vtk_celldata = build_vtk_grids(Val(format), mesh, node_set, n_visnodes, verbose,
                                                 output_directory, is_datafile, filename, Val(reinterpolate))

    # Interpolate data
    if is_datafile
      verbose && println("| Interpolating data...")
      if reinterpolate
        @timeit "interpolate data" interpolated_data = interpolate_data(Val(format),
                                                                        data, mesh,
                                                                        n_visnodes, verbose)
      else # Copy the raw solution data; only works for `vtu` format
        # Extract data shape information
        ndims_ = ndims(data) - 2
        n_variables = length(labels)
        # Save raw data as one 1D array for each variable
        @timeit "interpolate data" interpolated_data = reshape(data,
                                                               n_visnodes^ndims_ * n_elements,
                                                               n_variables)
      end
    end

    # Add data to file
    verbose && println("| Adding data to VTK file...")
    @timeit "add data to VTK file" begin
      if save_celldata
        add_celldata!(vtk_celldata, mesh, verbose)
      end

      # Only add data if it is a data file
      if is_datafile
        # Add solution variables
        for (variable_id, label) in enumerate(labels)
          verbose && println("| | Variable: $label...")
          @timeit label vtk_nodedata[label] = @views interpolated_data[:, variable_id]
        end

        if save_celldata
          # Add element variables
          for (label, variable) in element_variables
            verbose && println("| | Element variable: $label...")
            @timeit label vtk_celldata[label] = variable
          end

          # Add node variables
          for (label, variable) in node_variables
            verbose && println("| | Node variable: $label...")
            if reinterpolate
              # Show warning if node-level data of subcell limiting are reinterpolated.
              if label == "limiting_coefficient" && !has_warned_about_interpolation
                println("WARNING: The limiting coefficients are no continuous field but happens " *
                "to be represented by a piecewise-constant approximation. Thus, reinterpolation " *
                "does not give a meaningful representation.")
                has_warned_about_interpolation = true
              end
              @timeit "interpolate data" interpolated_cell_data = interpolate_data(Val(format),
                                                                    reshape(variable, size(variable)..., 1),
                                                                    mesh, n_visnodes, verbose)
            else
              @timeit "interpolate data" interpolated_cell_data = reshape(variable,
                                                                          n_visnodes^ndims_ * n_elements)
            end
            # Add to node_data
            @timeit label vtk_nodedata[label] = interpolated_cell_data
          end
        end
      end
    end

    # Save VTK file
    if is_datafile
      verbose && println("| Saving VTK file '$(vtk_nodedata.path)'...")
      @timeit "save VTK file" vtk_save(vtk_nodedata)
    end

    if save_celldata
      verbose && println("| Saving VTK file '$(vtk_celldata.path)'...")
      @timeit "save VTK file" vtk_save(vtk_celldata)
    end

    # Add to PVD file only if it is a datafile
    if !is_single_file
      if is_datafile
        verbose && println("| Adding to PVD file...")
        @timeit "add VTK to PVD file" begin
          pvd[time] = vtk_nodedata
          if save_celldata
            pvd_celldata[time] = vtk_celldata
          end
        end
        has_data = true
      else
        println("WARNING: file '$(filename)' will not be added to PVD file since it is a mesh file")
      end
    end

    # Update progress bar
    if !hide_progress
      next!(progress, showvalues=[(:finished, filename)])
    end
  end

  if !is_single_file
    # Save PVD file only if at least one data file was added
    if has_data
      verbose && println("| Saving PVD file '$(pvd_filename).pvd'...")
      @timeit "save PVD files" vtk_save(pvd)
    end

    if save_celldata
      verbose && println("| Saving PVD file '$(pvd_celldata_filename).pvd'...")
      @timeit "save PVD files" vtk_save(pvd_celldata)
    end
  end

  verbose && println("| done.\n")
  print_timer()
  println()
end


function assert_cells_elements(n_elements, mesh::TreeMesh, filename, meshfile)
  # Check if dimensions match
  if length(Trixi.leaf_cells(mesh.tree)) != n_elements
    error("number of elements in '$(filename)' do not match number of leaf cells in " *
          "'$(meshfile)' " *
          "(did you forget to clean your 'out/' directory between different runs?)")
  end
end


function assert_cells_elements(n_elements, mesh::StructuredMesh, filename, meshfile)
  # Check if dimensions match
  if prod(size(mesh)) != n_elements
    error("number of elements in '$(filename)' do not match number of cells in " *
          "'$(meshfile)' " *
          "(did you forget to clean your 'out/' directory between different runs?)")
  end
end


function assert_cells_elements(n_elements, mesh::UnstructuredMesh2D, filename, meshfile)
  # Check if dimensions match
  if length(mesh) != n_elements
    error("number of elements in '$(filename)' do not match number of cells in " *
          "'$(meshfile)' " *
          "(did you forget to clean your 'out/' directory between different runs?)")
  end
end


function assert_cells_elements(n_elements, mesh::Union{P4estMesh, T8codeMesh}, filename, meshfile)
  # Check if dimensions match
  if Trixi.ncells(mesh) != n_elements
    error("number of elements in '$(filename)' do not match number of cells in " *
          "'$(meshfile)' " *
          "(did you forget to clean your 'out/' directory between different runs?)")
  end
end


# default number of visualization nodes if a solution should be visualized
function get_default_nvisnodes_solution(nvisnodes, n_nodes, mesh::TreeMesh)
  if nvisnodes === nothing
    return 2 * n_nodes
  elseif nvisnodes == 0
    return n_nodes
  else
    return nvisnodes
  end
end

function get_default_nvisnodes_solution(nvisnodes, n_nodes,
                                        mesh::Union{StructuredMesh, UnstructuredMesh2D, P4estMesh, T8codeMesh})
  if nvisnodes === nothing || nvisnodes == 0
    return n_nodes
  else
    return nvisnodes
  end
end


# default number of visualization nodes if only the mesh should be visualized
function get_default_nvisnodes_mesh(nvisnodes, mesh::TreeMesh)
  if nvisnodes === nothing
    # for a Cartesian mesh, we do not need to interpolate
    return 1
  else
    return nvisnodes
  end
end

function get_default_nvisnodes_mesh(nvisnodes,
                                    mesh::Union{StructuredMesh, UnstructuredMesh2D, P4estMesh, T8codeMesh})
  if nvisnodes === nothing
    # for curved meshes, we need to get at least the vertices
    return 2
  else
    return nvisnodes
  end
end


function add_celldata!(vtk_celldata, mesh::TreeMesh, verbose)
  @timeit "add data to VTK file" begin
    leaf_cells = Trixi.leaf_cells(mesh.tree)
    # Add cell/element data to celldata VTK file
    verbose && println("| | cell_ids...")
    @timeit "cell_ids" vtk_celldata["cell_ids"] = leaf_cells
    verbose && println("| | element_ids...")
    @timeit "element_ids" vtk_celldata["element_ids"] = collect(1:length(leaf_cells))
    verbose && println("| | levels...")
    @timeit "levels" vtk_celldata["levels"] = mesh.tree.levels[leaf_cells]
  end

  return vtk_celldata
end


function add_celldata!(vtk_celldata, mesh::StructuredMesh, verbose)
  @timeit "add data to VTK file" begin
    # Add element data to celldata VTK file
    verbose && println("| | element_ids...")
    @timeit "element_ids" vtk_celldata["element_ids"] = collect(1:prod(mesh.cells_per_dimension))
  end

  return vtk_celldata
end


function add_celldata!(vtk_celldata, mesh::UnstructuredMesh2D, verbose)
  @timeit "add data to VTK file" begin
    # Add element data to celldata VTK file
    verbose && println("| | element_ids...")
    @timeit "element_ids" vtk_celldata["element_ids"] = collect(1:length(mesh))
  end

  return vtk_celldata
end


function add_celldata!(vtk_celldata, mesh::P4estMesh, verbose)
  # Create temporary storage for the tree_ids and levels.
  tree_ids = zeros( Trixi.ncells(mesh) )
  cell_levels = zeros( Trixi.ncells(mesh) )
  # Set global counters.
  tree_counter = 1
  cell_counter = 1
  # Iterate through the p4est trees and each of their quadrants.
  # Assigns the tree index values. Also, grab and assign the level value.
  for tree in Trixi.unsafe_wrap_sc(Trixi.P4est.p4est_tree_t, unsafe_load(mesh.p4est).trees)
    for quadrant in Trixi.unsafe_wrap_sc(Trixi.P4est.p4est_quadrant_t, tree.quadrants)
      tree_ids[cell_counter] = tree_counter
      cell_levels[cell_counter] = quadrant.level
      cell_counter += 1
    end
    tree_counter += 1
  end
  @timeit "add data to VTK file" begin
    # Add tree/element data to celldata VTK file
    verbose && println("| | tree_ids...")
    @timeit "tree_ids" vtk_celldata["tree_ids"] = tree_ids
    verbose && println("| | element_ids...")
    @timeit "element_ids" vtk_celldata["element_ids"] = collect(1:Trixi.ncells(mesh))
    verbose && println("| | levels...")
    @timeit "levels" vtk_celldata["levels"] = cell_levels
  end

  return vtk_celldata
end

function add_celldata!(vtk_celldata, mesh::T8codeMesh, verbose)
  # Create temporary storage for the tree_ids and levels.
  tree_ids = zeros( Trixi.ncells(mesh) )

  elem_counter = 1
  num_local_trees = Trixi.t8_forest_get_num_local_trees(mesh.forest)
  for itree in 1:num_local_trees
      num_elements_in_tree = Trixi.t8_forest_get_tree_num_elements(mesh.forest, itree-1)
      for ielement in 1:num_elements_in_tree
          tree_ids[elem_counter] = itree
          elem_counter += 1
      end
  end

  levels = Trixi.trixi_t8_get_local_element_levels(mesh.forest)

  @timeit "add data to VTK file" begin
    # Add tree/element data to celldata VTK file
    verbose && println("| | tree_ids...")
    @timeit "tree_ids" vtk_celldata["tree_ids"] = tree_ids
    verbose && println("| | element_ids...")
    @timeit "element_ids" vtk_celldata["element_ids"] = collect(1:Trixi.ncells(mesh))
    verbose && println("| | levels...")
    @timeit "levels" vtk_celldata["levels"] = levels
  end

  return vtk_celldata
end

function expand_filename_patterns(patterns)
  filenames = String[]

  for pattern in patterns
    if startswith(pattern, '/')
      # Glob.glob does not support absolute paths
      append!(filenames, glob(relpath(pattern)))
    else
      append!(filenames, glob(pattern))
    end
  end

  return filenames
end

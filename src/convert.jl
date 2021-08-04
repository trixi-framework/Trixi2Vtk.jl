"""
    trixi2vtk(filename::AbstractString...;
              format=:vtu, verbose=false, hide_progress=false, pvd=nothing,
              output_directory=".", nvisnodes=nothing)

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
                   (the default `nothing` is converted to `false`
                   for `StructuredMesh`/`UnstructuredMesh2D` and `true` for `TreeMesh`)

# Examples
```julia
julia> trixi2vtk("out/solution_000*.h5")
[...]
```
"""
function trixi2vtk(filename::AbstractString...;
                   format=:vtu, verbose=false, hide_progress=false, pvd=nothing,
                   output_directory=".", nvisnodes=nothing, save_celldata=nothing)
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
    progress = Progress(length(filenames), 0.5, "Converting .h5 to .$(format)...", 40)
  end

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

    if save_celldata === nothing
      # If no value for `save_celldata` is specified,
      # use true for TreeMesh and false for StructuredMesh or UnstructuredMesh2D
      save_celldata = isa(mesh, TreeMesh)
    end

    # Read data only if it is a data file
    if is_datafile
      verbose && println("| Reading data file...")
      @timeit "read data" (labels, data, n_elements, n_nodes,
                           element_variables, time) = read_datafile(filename)

      assert_cells_elements(n_elements, mesh, filename, meshfile)

      # Determine resolution for data interpolation
      n_visnodes = get_default_nvisnodes(nvisnodes, n_nodes, mesh)
    else
      # If file is a mesh file, do not interpolate data
      n_visnodes = 1
    end

    # Create output directory if it does not exist
    mkpath(output_directory)

    # Build VTK grids
    vtk_nodedata, vtk_celldata = build_vtk_grids(Val(format), mesh, n_visnodes, verbose,
                                                 output_directory, is_datafile, filename)

    # Interpolate data
    if is_datafile
      verbose && println("| Interpolating data...")
      @timeit "interpolate data" interpolated_data = interpolate_data(Val(format),
                                                                      data, mesh,
                                                                      n_visnodes, verbose)
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


function assert_cells_elements(n_elements, mesh::P4estMesh, filename, meshfile)
  # Check if dimensions match
  if Trixi.ncells(mesh) != n_elements
    error("number of elements in '$(filename)' do not match number of cells in " *
          "'$(meshfile)' " *
          "(did you forget to clean your 'out/' directory between different runs?)")
  end
end


function get_default_nvisnodes(nvisnodes, n_nodes, mesh::TreeMesh)
  if nvisnodes === nothing
    return 2 * n_nodes
  elseif nvisnodes == 0
    return n_nodes
  else
    return nvisnodes
  end
end


function get_default_nvisnodes(nvisnodes, n_nodes,
                               mesh::Union{StructuredMesh, UnstructuredMesh2D, P4estMesh})
  if nvisnodes === nothing || nvisnodes == 0
    return n_nodes
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

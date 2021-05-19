
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


# Calculate and interpolation matrix (Vandermonde matrix) between two given sets of nodes
function polynomial_interpolation_matrix(nodes_in, nodes_out)
  n_nodes_in = length(nodes_in)
  n_nodes_out = length(nodes_out)
  wbary_in = barycentric_weights(nodes_in)
  vdm = zeros(n_nodes_out, n_nodes_in)

  for k = 1:n_nodes_out
    match = false
    for j = 1:n_nodes_in
      if isapprox(nodes_out[k], nodes_in[j], rtol=eps())
        match = true
        vdm[k, j] = 1
      end
    end

    if match == false
      s = 0.0
      for j = 1:n_nodes_in
        t = wbary_in[j] / (nodes_out[k] - nodes_in[j])
        vdm[k, j] = t
        s += t
      end
      for j = 1:n_nodes_in
        vdm[k, j] = vdm[k, j] / s
      end
    end
  end

  return vdm
end


# Calculate the barycentric weights for a given node distribution.
function barycentric_weights(nodes)
  n_nodes = length(nodes)
  weights = ones(n_nodes)

  for j = 2:n_nodes, k = 1:(j-1)
    weights[k] *= nodes[k] - nodes[j]
    weights[j] *= nodes[j] - nodes[k]
  end

  for j = 1:n_nodes
    weights[j] = 1 / weights[j]
  end

  return weights
end


# From FLUXO (but really from blue book by Kopriva)
function gauss_lobatto_nodes_weights(n_nodes::Integer)
  # From Kopriva's book
  n_iterations = 10
  tolerance = 1e-15

  # Initialize output
  nodes = zeros(n_nodes)
  weights = zeros(n_nodes)

  # Get polynomial degree for convenience
  N = n_nodes - 1

  # Calculate values at boundary
  nodes[1] = -1.0
  nodes[end] = 1.0
  weights[1] = 2 / (N * (N + 1))
  weights[end] = weights[1]

  # Calculate interior values
  if N > 1
    cont1 = pi/N
    cont2 = 3/(8 * N * pi)

    # Use symmetry -> only left side is computed
    for i in 1:(div(N + 1, 2) - 1)
      # Calculate node
      # Initial guess for Newton method
      nodes[i+1] = -cos(cont1*(i+0.25) - cont2/(i+0.25))

      # Newton iteration to find root of Legendre polynomial (= integration node)
      for k in 0:n_iterations
        q, qder, _ = calc_q_and_l(N, nodes[i+1])
        dx = -q/qder
        nodes[i+1] += dx
        if abs(dx) < tolerance * abs(nodes[i+1])
          break
        end
      end

      # Calculate weight
      _, _, L = calc_q_and_l(N, nodes[i+1])
      weights[i+1] = weights[1] / L^2

      # Set nodes and weights according to symmetry properties
      nodes[N+1-i] = -nodes[i+1]
      weights[N+1-i] = weights[i+1]
    end
  end

  # If odd number of nodes, set center node to origin (= 0.0) and calculate weight
  if n_nodes % 2 == 1
    _, _, L = calc_q_and_l(N, 0)
    nodes[div(N, 2) + 1] = 0.0
    weights[div(N, 2) + 1] = weights[1] / L^2
  end

  return nodes, weights
end


# From FLUXO (but really from blue book by Kopriva)
function calc_q_and_l(N::Integer, x::Float64)
  L_Nm2 = 1.0
  L_Nm1 = x
  Lder_Nm2 = 0.0
  Lder_Nm1 = 1.0

  local L
  for i in 2:N
    L = ((2 * i - 1) * x * L_Nm1 - (i - 1) * L_Nm2) / i
    Lder = Lder_Nm2 + (2 * i - 1) * L_Nm1
    L_Nm2 = L_Nm1
    L_Nm1 = L
    Lder_Nm2 = Lder_Nm1
    Lder_Nm1 = Lder
  end

  q = (2 * N + 1)/(N + 1) * (x * L - L_Nm2)
  qder = (2 * N + 1) * L

  return q, qder, L
end
calc_q_and_l(N::Integer, x::Real) = calc_q_and_l(N, convert(Float64, x))

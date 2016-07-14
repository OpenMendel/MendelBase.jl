################################################################################
# This set of functions houses the Search optimization algorithm.
################################################################################

# using DataStructures
 
export set_parameter_defaults, optimize

"""
This function constructs the parameter structure.
"""
function set_parameter_defaults(keyword::Dict{ASCIIString, Any})
  #
  cases = keyword["cases"]
  constraints = keyword["constraints"]
  goal = keyword["goal"]
  output_unit = keyword["output_unit"]
  pars = keyword["parameters"]
  points = keyword["points"]
  standard_errors = keyword["standard_errors"]
  travel = keyword["travel"]
  title = keyword["title"]
  #
  # Check that input variables are properly defined.
  #
  if constraints < 0
    throw(ArgumentError(
      "The number of constraints ($constraints) must be nonnegative.\n \n"))
  end
  if pars < 0
    throw(ArgumentError(
      "The number of parameters ($pars) must be nonnegative.\n \n"))
  end
  if points < 0
    throw(ArgumentError(
      "The number of grid points ($points) must be nonnegative.\n \n"))
  end
  if goal != "minimize" && goal != "maximize"
    throw(ArgumentError(
      "The variable goal ($goal) must be minimize or maximize.\n \n"))
  end
  if travel != "grid" && travel != "search"
    throw(ArgumentError(
      "The variable travel ($travel) must be grid or search.\n \n"))
  end
  #
  # Name the parameters par 1, par 2, and so forth.
  #
  pname = ["par" for i = 1:pars]
  for i = 1:pars
    pname[i] = pname[i]*" $i"
    pname[i] = rpad(pname[i], 8, ' ')
  end
  #
  # Set other array defaults.
  #
  grid = zeros(points, pars)
  par = zeros(pars)
  pmin = zeros(pars)
  pmin[1:pars] = -Inf
  pmax = zeros(pars)
  pmax[1:pars] = Inf
  constraint = zeros(constraints, pars)
  constraint_level = zeros(constraints)
  level = zeros(constraints)
  if travel == "grid"
    function_value = zeros(points)
  else
    function_value = zeros(2)
  end
  parameter = Parameter(cases, constraints, goal, output_unit, pars, 
    points, standard_errors, title, travel, constraint, constraint_level, 
    function_value, grid, pmin, pmax, pname, par)
  #
  # Let the user change the defaults.
  #
  return parameter
end # function set_parameter_defaults

"""
This function computes the values of a function f over a user
defined grid of points or the constrained optimum of f by the
method of recursive quadratic programming. The hessian d2f of
f is approximated by Davidon's secant update unless an exact
value is supplied by the user. Users must supply the objective
function fun. Parameter initial values and constraints are
passed through the parameter data structure.
"""
function optimize(fun::Function, parameter::Parameter)
  #
  # Set variables controlling convergence and numerical differentiation.
  #
  conv_crit = 1e-6
  conv_tests = 4
  max_iters = 1000
  max_steps = 3
  max_step_length = Inf
  dp = 1e-7
  small = 1e-5
  #
  # Write a header.
  #
  cases = parameter.cases
  io = parameter.output_unit
  travel = parameter.travel
  goal = parameter.goal
  println(io, " \n                Search, Julia Version\n")
  println(io, "            (c) Copyright Kenneth Lange, 2016\n")
  if parameter.title != ""
    println(io, "Title: " * parameter.title)
  end
  println(io, "Grid or search option: " * travel)
  println(io, "Minimize or maximize: " * goal)
  #
  # Abbreviate the names of passed variables.
  #
  constraints = parameter.constraints
  pars = parameter.parameters
  points = parameter.points
  travel = parameter.travel
  pname = parameter.name
  par = parameter.par
  pmin = parameter.min
  pmax = parameter.max
  constraint = parameter.constraint
  level = parameter.constraint_level
  function_value = parameter.function_value
  #
  # Create variables required by searching.
  #
  grad_provided = false
  hess_provided = false
  f = 0.0
  df = zeros(pars)
  grad = zeros(pars)
  d2f = zeros(pars, pars)
  hess = zeros(pars, pars)
  df_old = zeros(pars)
  par_old = zeros(pars)
  #
  # Check for bound and constraint errors.
  #
  error = check_constraints(constraint, level, pname, par, pmin, pmax, io, travel)
  if error; return (par, 0.0); end
  #
  # Compute function values over a user defined grid of points.
  # Keep track of the optimal point and optimal function value.
  #
  fmin = Inf
  iter_min = 0
  best_point = zeros(pars)
  if travel == "grid"
    for iteration = 1:points
      parameter.par = vec(parameter.grid[iteration, :])
      (f, df, d2f) = fun(parameter.par)
      iteration_output!(io, f, iteration, parameter.par, pname)
      function_value[iteration] = f
      if goal == "maximize"
        f = - f
      end
      if f < fmin
        fmin = f
        iter_min = iteration
        copy!(best_point, par)
      end
    end
    #
    # Return the optimal point and optimal function value.
    #
    if goal == "minimize"
      println(io, " \nThe minimum function value of ", round(fmin, 5),
        " occurs at iteration ", iter_min, ".")
      return (best_point, fmin)
    else
      println(io, " \nThe maximum function value of ", round(- fmin, 5),
        " occurs at iteration ", iter_min, ".")
      return (best_point, - fmin)
    end
    #
    # Otherwise minimize the objective function. df and d2f
    # are the first two differentials. ker is the kernel of
    # the constraint matrix.
    #
  else
    df = zeros(pars)
    d2f = eye(pars, pars)
    (matrix_u, d, matrix_v) = svd(constraint, thin = false)
    ker = matrix_v[:, constraints + 1:pars]
    conv_count = 0
    forward = true
    #
    # Fetch the first function value and the two differentials.
    # Output the first iteration.
    #
    iteration = 1
    (f, grad_provided, grad, hess) = dfun(fun, parameter, dp, forward)
    iteration_output!(io, f, iteration, par, pname)
    function_value[1] = f
    #
    # Copy the user supplied first and second differentials.
    #
    hess_provided =  typeof(hess) != Void
    if goal == "minimize"
      copy!(df, grad)
      if hess_provided; copy!(d2f, hess); end
    else
      f = - f
      copy!(df, - grad)
      if hess_provided; copy!(d2f, - hess); end    
    end    
    #
    # Preserve the current point.
    #
    copy!(best_point, par)
    fmin = f
    iter_min = 1
    #
    # Enter the main iteration loop.
    #
    for iteration = 2:max_iters
    #
    # Solve the quadratic programming problem. If cycles = 0,
    # then the quadratic programming algorithm has failed. Reset
    # the Hessian approximation to the identity and try again.
    #
      for j = 1:2
        tableau = create_tableau(d2f, df, constraint, level, par)
        (cycles, delta) =
        quadratic_program(tableau, par, pmin, pmax, pars, constraints)
        #
        # Check whether quadratic programming succeeded.
        #
        if cycles > 0
          break
        else
          c = Inf
          for i = 1:pars
            if d2f[i, i] > 0.0; c = min(c, d2f[i, i]); end
          end
          #
          # In an emergency revert to steepest descent.
          #
          if c >= Inf
            d2f = eye(pars, pars)
          #
          # Otherwise, bump up the diagonal entries of the hessian.
          #
          else
            for i = 1:pars
              if d2f[i, i] <= c; d2f[i, i] = c * (1.0 + rand()); end
            end
          end
        end
      end
      #
      # Compute a new function value. If it represents a sufficient
      # decrease, then accept the new point. Otherwise backtrack by
      # fitting a quadratic through the two current points with the
      # given slope at the first point. Quit, if necessary, after the
      # maximum number of steps.
      #
      l2_norm = norm(delta)
      if l2_norm > max_step_length
        delta = max_step_length * delta / l2_norm
      end
      d = min(dot(df, delta), 0.0)
      copy!(df_old, df)
      copy!(par_old, par)
      f_old = f
      t = 1.0
      steps = -1
      for step = 0:max_steps
        steps = steps + 1
        #
        # Fetch another function value, and if provided, the user
        # supplied second differential.
        #
        par = par_old + t * delta
        parameter.par = par
        (f, grad_provided, grad, hess) = dfun(fun, parameter, dp, forward)
        if goal == "minimize"
          copy!(df, grad)
          if hess_provided; copy!(d2f, hess); end
        else
          f = - f
          copy!(df, - grad)
          if hess_provided; copy!(d2f, - hess); end    
        end
        #
        # Test for a sufficient decrease in the objective.
        #
        if f <= f_old + t * d / 10.0 || step == max_steps
          break
        else
          t = max( -d * t * t / (2.0 * (f - f_old - t * d)), t / 10.0)
        end
      end
      #
      # Output the current iteration, and update the best point.
      #
      if goal == "minimize"
        iteration_output!(io, f, iteration, par, pname, steps)
      else
        iteration_output!(io, - f, iteration, par, pname, steps)
      end
      if f < fmin
        fmin = f
        iter_min = iteration      
        copy!(best_point, par)
      end
      #
      # Check the convergence criterion. If it has been satisfied a
      # sufficient number of times, then exit the main loop. Otherwise,
      # output the current iteration and continue the search.
      #
      if abs(f_old - f) > conv_crit
        forward = true
        conv_count = 0
      else
        forward = false
        conv_count = conv_count + 1
      end
      if conv_count >= conv_tests
        break
      end
      #
      # Update the approximate Hessian by Davidon's secant method.
      #
      if !hess_provided
        delta = par - par_old
        u = df - df_old - d2f * delta
        c = dot(u, delta)
        if c^2 > (1e-8) * dot(u, u) * dot(delta, delta)
          c = -1.0 / c
          #
          # Ensure that the updated Hessian is positive definite on
          # the kernel of the constraint matrix.
          #
          matrix = inv(ker' * d2f * ker)
          v = ker' * u
          d = dot(v, matrix * v)
          c = min(c, 0.8 / d)
          d2f = d2f - c * u * u'
        end
      end
    end
  end
  #
  # Output the optimal function value.
  #
  if goal == "minimize"
    println(io, " \nThe minimum function value of ", round(fmin, 5),
      " occurs at iteration ", iter_min, ".")
    function_value[2] = fmin
  else
    println(io, " \nThe maximum function value of ", round(- fmin, 5),
      " occurs at iteration ", iter_min, ".")
    function_value[2] = - fmin
  end
  #
  # If the asymptotic covariance matrix is desired, then
  # record which parameters occur on a boundary
  #
  boundary = falses(pars)
  if parameter.standard_errors
    for i = 1:pars
      boundary[i] = (par[i] <= pmin[i] + small || par[i] >= pmax[i] - small)
    end
    boundaries = sum(boundary)
    #
    # If the exact first differential is available, approximate the
    # columns of the second differential by forward differences of
    # the first differential.
    #
    if !hess_provided
      if grad_provided
        for i = 1:pars
          if boundary[i]; continue; end
          d = dp * max(abs(par[i]), 1.0)
          if par[i] + d >= pmax[i]; d = -d; end
          parameter.par[i] = par[i] + d
          (g, grad_provided, grad, hess) = dfun(fun, parameter, dp, forward)
          if goal == "maximize"; g = - g; end
          parameter.par[i] = parameter.par[i] - d
          d2f[1:i, i] = (grad[1:i] - df[1:i]) / d
          d2f[i, 1:i] = d2f[1:i, i]
        end
      #
      # Otherwise, use the second difference formula of Dennis and Schnabel
      # to approximate the entries of the second differential. See page
      # on page 321 of their book.
      #
      else
        #
        # If necessary, recompute the first differential using central differences.
        #
        if forward
          forward = false
          parameter.par = par
          (f, grad_provided, df, hess) = dfun(fun, parameter, dp, forward)
        end
        #
        # Compute and store function values and parameter increments.
        #
        u = zeros(pars)
        delta = zeros(pars)
        dp23 = cbrt(dp)^2
        for i = 1:pars
          if boundary[i]; continue; end
          d = dp23 * max(abs(par[i]), 1.0)
          if parameter.par[i] + d >= pmax[i]; d = -d; end
          parameter.par[i] = parameter.par[i] + d
          (g, grad, hess) = fun(parameter.par)
          if goal == "maximize"; g = - g; end
          parameter.par[i] = parameter.par[i] - d
          delta[i] = d
          u[i] = g
        end
        #
        # Apply these to compute the entries of the second differential.
        #
        for i = 1:pars
          for j = 1:i
            if !boundary[i] && !boundary[j]
              parameter.par[i] = parameter.par[i] + delta[i]
              parameter.par[j] = parameter.par[j] + delta[j]
              (g, grad, hess) = fun(parameter.par)
              if goal == "maximize"; g = - g; end
              parameter.par[i] = parameter.par[i] - delta[i]
              parameter.par[j] = parameter.par[j] - delta[j]
              d2f[i, j] = (g - u[i] - u[j] + f) / (delta[i] * delta[j])
              d2f[j, i] = d2f[i, j]
            end
          end
        end
      end
    end
    #
    # Adjust the Hessian for a least squares problem by dividing
    # by the residual mean square.
    #
    reduced = pars - constraints - boundaries
    if cases > 0
      n = cases - reduced
      if n > 0 && f > 0.0
        sigma_sq = 2.0 * f / n
        d2f = d2f / sigma_sq
      else
        fill!(d2f, 0.0)
      end
    end
    #
    # After these preliminaries, compute the asymptotic covariances.
    #
    if reduced > 0
      asy_cov = asymptotic_covariances(io, constraint, boundary, d2f, pname)
    else
      write(io, " \n The asymptotic covariance matrix is undefined.")
      println(io, "")
    end
  end
  #
  # Return the optimal point and optimal function value.
  #
  if goal == "minimize"
    return (best_point, fmin)
  else
    return (best_point, - fmin)
  end
end # function optimize

"""
This function calculates the tableau used in minimizing
the quadratic 0.5 x' Q x + r' x, subject to Ax = b and
parameter lower and upper bounds.
"""
function create_tableau(matrix_q::Matrix{Float64}, r::Vector{Float64},
  matrix_a::Matrix{Float64}, b::Vector{Float64}, x::Vector{Float64})

  m = size(matrix_a, 1)
  #
  # Create the tableau in the absence of constraints.
  #
  if m == 0
    tableau = [matrix_q (-r); -r' 0]
  else
  #
  # In the presence of constraints compute a constant mu via
  # the Gerschgorin circle theorem so that Q + mu * A' * A
  # is positive definite.
  #
    (matrix_u, d, matrix_v) = svd(matrix_a, thin = false)
    matrix_p = matrix_v * matrix_q * matrix_v'
    mu = 0.0
    for i = 1:m
      mu = max((norm(matrix_p[:, i], 1) - 2.0 * matrix_p[i, i]) / d[i]^2, mu)
    end
    mu = 2.0 * mu
    #
    # Now create the tableau.
    #
    tableau = [matrix_q + mu * matrix_a' * matrix_a matrix_a' (-r);
               matrix_a zeros(m, m) b - matrix_a * x;
               -r' (b - matrix_a * x)' 0]
  end
  return tableau
end # function create_tableau

"""
This subroutine solves the p-dimensional quadratic programming problem

 min [df * delta + 0.5 * delta' * d^2 f * delta]
 subject to: constraint * delta = 0 and pmin <= par + delta <= pmax.

See: Jennrich JI, Sampson PF (1978) Some problems faced in making
a variance component algorithm into a general mixed model program.
Proceedings of the Eleventh Annual Symposium on the Interface.
Gallant AR, Gerig TM, editors. Institute of Statistics,
North Carolina State University.
"""
function quadratic_program(tableau::Matrix{Float64}, par::Vector{Float64},
  pmin::Vector{Float64}, pmax::Vector{Float64}, p::Int, c::Int)

  delta = zeros(par)
  #
  # See function create_tableau for the construction of the tableau.
  # For checking tolerance, set diag to the diagonal elements of tableau.
  # Begin by sweeping on those diagonal elements of tableau corresponding
  # to the parameters. Then sweep on the diagonal elements corresponding
  # to the constraints. If any parameter fails the tolerance test, then
  # return and reset the approximate Hessian.
  #
  small = 1e-5
  tol = 1e-8
  d = diag(tableau)
  for i = 1:p
    if d[i] <= 0.0 || tableau[i, i] < d[i] * tol
      return (0, delta)
    else
      sweep!(tableau, i, false)
    end
  end
  swept = trues(p)
  for i = p + 1:p + c
    if tableau[i, i] >= 0.0
      return (0, delta)
    else
      sweep!(tableau, i, false)
    end
  end
  #
  # Take a step in the direction tableau(i, end) for the parameters i
  # that are currently swept. If a boundary is encountered, determine
  # the maximal fractional step possible.
  #
  cycle_main_loop = false
  for iteration = 1:1000
    a = 1.0
    for i = 1:p
      if swept[i]
        ui = tableau[i, end]
        if ui > 0.0
          ai = pmax[i] - par[i] - delta[i]
        else
          ai = pmin[i] - par[i] - delta[i]
        end
        if abs(ui) > 1e-10
          a = min(a, ai / ui)
        end
      end
    end
    #
    # Take the fractional step for the currently swept parameters, and
    # reset the transformed partial derivatives for these parameters.
    #
    for i = 1:p
      if swept[i]
        ui = tableau[i, end]
        delta[i] = delta[i] + a * ui
        tableau[i, end] = (1.0 - a) * ui
        tableau[end, i] = tableau[i, end]
      end
    end
    #
    # Find a swept parameter that is critical, and inverse sweep it.
    # Go back and try to take another step or fractional step.
    #
    cycle_main_loop = false
    for i = 1:p
      critical = pmin[i] >= par[i] + delta[i] - small
      critical = critical || pmax[i]<= par[i] + delta[i] + small
      if swept[i] && abs(tableau[i, i])>1e-10 && critical
        sweep!(tableau, i, true)
        swept[i] = false
        cycle_main_loop = true
        break
      end
    end
    if cycle_main_loop; continue; end
    #
    # Find an unswept parameter that violates the KKT condition
    # and sweep it. Go back and try to take a step or fractional step.
    # If no such parameter exists, then the problem is solved.
    #
    for i = 1:p
      ui = tableau[i, end]
      violation = ui > 0.0 && pmin[i] >= par[i] + delta[i] - small
      violation = violation || (ui<0.0 && pmax[i]<= par[i] + delta[i] + small)
      if !swept[i] && violation
        sweep!(tableau, i, false)
        swept[i] = true
        cycle_main_loop = true
        break
      end
    end
    if cycle_main_loop; continue; end
    return (iteration, delta)
  end
  return (0, delta)
end # function quadratic_program

"""
This function sweeps or inverse sweeps the symmetric tableau A
on its kth diagonal entry.
"""
function sweep!(matrix_a::Matrix{Float64}, k::Int, inverse::Bool = false)

  p = 1.0 / matrix_a[k, k]
  v = matrix_a[:, k]
  matrix_a[:, k] = 0.0
  matrix_a[k, :] = 0.0
  if inverse
    v[k] = 1.0
  else
    v[k] = -1.0
  end
  for i = 1:size(matrix_a, 1)
    pv = p * v[i]
    matrix_a[:, i] = matrix_a[:, i] - pv * v
  end
end # function sweep!

"""
This function controls the computation of f and its first two differentials.
When the first partials are computed numerically, forward or central 
differences are used depending on the logical variable forward. The numerical
differentiation interval is adjusted to take into account the magnitude
of a parameter and whether it lies on its upper bound.
"""
function dfun(fun::Function, parameter::Parameter, dp::Float64, forward::Bool)

  par = parameter.par
  pmax = parameter.max
  (f, df, d2f) = fun(par)
  grad_provided = typeof(df) != Void
  if !grad_provided
    df = zeros(length(par))
    for i = 1:length(par)
      d = dp * max(abs(par[i]), 1.0)
      psave = par[i]
      if forward
        if psave + d > pmax[i]; d = -d; end
        par[i] = par[i] + d
        (fplus, grad, hess) = fun(par)
        df[i] = (fplus - f) / d
      else
        par[i] = par[i] + d
        (fplus, grad, hess) = fun(par)
        par[i] = psave - d
        (fminus, grad, hess) = fun(par)
        df[i] = (fplus - fminus) / (d + d)
      end
      par[i] = psave
    end
  end
  return (f, grad_provided, df, d2f)
end # function dfun

"""
This function checks for constraint violations.
The returned value conveys the nature of the violation.
"""
function check_constraints(constraint::Matrix{Float64},
  level::Vector{Float64}, pname::Vector{ASCIIString}, par::Vector{Float64},
  pmin::Vector{Float64}, pmax::Vector{Float64}, io::IO, travel::ASCIIString)

  tol = 1e-4
  error = false
  (constraints, pars) = size(constraint)
  #
  # Check for bound violations.
  #
  for i = 1:pars
    if par[i] < pmin[i]
      println(io,
        " \nError: Parameter ", i, " is less than its minimum.")
      error = true
    elseif par[i] > pmax[i]
      println(io,
        " \nError: Parameter ", i, " is greater than its maximum.")
      error = true
    end
    if pmin[i] > pmax[i] - tol
      println(io, " \nError: The bounds on parameter ", i,
                  " are too close or inconsistent.")
      error = true
    end
  end
  #
  # Check for constraint violations.
  #
  if constraints > 0
    y = constraint * par
    for i = 1:constraints
      if abs(y[i] - level[i]) > tol
        println(io,
          " \nError: Linear equality constraint ", i, " is not satisfied.")
        error = true
      end
      if countnz(constraint[i, :]) != 1
        continue
      else
        for j = 1:pars
          if constraint[i, j] != 0.0 &&
          (abs(par[j] - pmin[j]) < tol || abs(par[j] - pmax[j]) < tol)
            println(io,
            " \nError: parameter ", j, " is constrained to one of its bounds.")
            error = true
          end
        end
      end
    end
    if rank(constraint) < size(constraint, 1)
      println(io, " \nError: Some equality constraints are redundant.")
      error = true
    end
  end
  #
  # Echo the bounds and the constraints.
  #
  if travel == "search"
    println(io, " \nParameter minima and maxima:")
    name = join(pname, "    ")
    println(io, " \n    ", name)
    @printf(io, " \n%12.4e", pmin[1])
    for i = 2:pars
      @printf(io, "%12.4e", pmin[i])
    end
    @printf(io, " \n%12.4e", pmax[1])
    for i = 2:pars
      @printf(io, "%12.4e", pmax[i])
    end
    println(io, " ")
    if constraints > 0
      println(io, " \nParameter constraints:")
      println(io, " \n    ", name, "    level ")
      println(io, "")
      for i = 1:constraints
        for j = 1:pars
          @printf(io, "%12.4e", constraint[i, j])
        end
        @printf(io, "%12.4e\n", level[i])
      end
    end
  end
  return error
end # function check_constraints

"""
This function outputs the current iteration.
"""
function iteration_output!(io::IO, f::Float64, iteration::Int,
  par::Vector{Float64}, pname::Vector{ASCIIString}, steps::Int = 0)

  if iteration == 1
    name = join(pname, "    ")
    println(io, " \niter  steps   function    ", name)
  end
  @printf(io, " \n %2i   %2i   %11.4e ", iteration, steps, f)
  for i = 1:length(par)
     @printf(io, "%12.4e", par[i])
  end
  println(io, " ")
end # function iteration_output!

"""
This function computes the asymptotic covariance matrix of the
parameter estimates.
"""
function asymptotic_covariances(io::IO, constraint::Matrix{Float64},
  boundary::BitArray{1}, d2f::Matrix{Float64}, pname::Vector{ASCIIString})
  #
  # Reparameterize to find the asymptotic covariance matrix.
  # Add an extra constraint for each parameter occurring on a boundary.
  #
  (constraints, pars) = size(constraint)
  asy_cov = zeros(pars, pars)
  asy_cov[1:constraints, :] = constraint
  j = constraints + 1
  for i = 1:pars
    if boundary[i]
      asy_cov[j, i] = 0.0
      j = j + 1
    end
  end
  (matrix_u, d, matrix_v) = svd(asy_cov, thin = false)
  ker = matrix_v[:, j:pars]
  matrix = ker' * d2f * ker
  if isposdef(matrix)
    asy_cov = ker * inv(matrix) * ker'
  else
    println(io, " \nThe asymptotic covariance matrix cannot be computed.\n")
    println(io, " ")
    return nothing
  end
  #
  # Standardize the asymptotic covariance matrix.
  #
  for i = 1:pars
    if asy_cov[i, i] <= 1e-10
      asy_cov[i, i] = 0.0
    else
      asy_cov[i, i] = sqrt(asy_cov[i, i])
    end
    for j = 1:i - 1
      if asy_cov[i, i] > 1e-5 &&  asy_cov[j, j] > 1e-5
        asy_cov[i, j] = asy_cov[i, j] / (asy_cov[i, i] * asy_cov[j, j])
      else
        asy_cov[i, j] = 0.0
        asy_cov[j, i] = 0.0
      end
    end
  end
  #
  # Output the results.
  #
  println(io, " \nThe asymptotic standard errors of the parameters:")
  name = join(pname, "    ")
  println(io, " \n    ", name)
  println(io, " ")
  for i = 1:pars
    @printf(io, "%12.4e", asy_cov[i, i])
  end
  println(io, "")
  println(io, " \nThe asymptotic correlation matrix of the parameters:")
  println(io, " \n    ", name)
  println(io, "")
  @printf(io, "%10.4f\n", 1.0)
  for i = 2:pars
    println(io, "")
    @printf(io, "%10.4f", asy_cov[i, 1])
    for j = 2:i - 1
      @printf(io, "%12.4f", asy_cov[i, j])
    end
    if asy_cov[i, i] > 1e-5
      @printf(io, "%12.4f", 1.0)
    else
      @printf(io, "%12.4f", 0.0)
    end
    println(io, "")
  end
end # function asymptotic_covariances


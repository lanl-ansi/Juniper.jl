# Juniper

Status:
[![CI](https://github.com/lanl-ansi/Juniper.jl/workflows/CI/badge.svg)](https://github.com/lanl-ansi/Juniper.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/lanl-ansi/Juniper.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/Juniper.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://lanl-ansi.github.io/Juniper.jl/stable)
</p>

# The Idea

You have a nonlinear problem with discrete variables (MINLP) and want some more control over the branch and bound part.
The relaxation should be solveable by any solver you prefer. Some solvers might not be able to solve the mixed integer part by themselves.

Juniper (Jump Nonlinear Integer Program solver) is a heuristic for optimization problems with non-convex functions.
If you need the global optimum, check out [Alpine](http://github.com/lanl-ansi/Alpine.jl).

# Basic usage

Juniper can be installed via the Julia package manager,

`] add JuMP, Juniper`

Add it to your project with,

```julia
using JuMP, Juniper
```

You will also have to have an NLP solver for setting up Juniper (e.g., [Ipopt](https://github.com/jump-dev/Ipopt.jl)), 

```julia
using Ipopt
```

Define a Juniper optimizer,

```julia
nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver)
```
The provided `nl_solver` is used by Juniper to solve continuous nonlinear sub-problems while Juniper searches for acceptable assignments to the discrete variables.

Give Juniper a try:

```julia
import LinearAlgebra: dot
m = Model(minlp_solver)

v = [10,20,12,23,42]
w = [12,45,12,22,21]
@variable(m, x[1:5], Bin)

@objective(m, Max, dot(v,x))

@constraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)

optimize!(m)

# retrieve the objective value, corresponding x values and the solver status
println(termination_status(m))
println(objective_value(m))
println(value.(x))
```

To solve problems with more complex nonlinear functions use the `@NLconstraint` and `@NLobjective` features of JuMP models.

If Juniper has difficulty finding feasible solutions on your model, try adding a mip solver (e.g., [HiGHS](https://github.com/jump-dev/HiGHS.jl)) to run a _feasiblity pump_,

```julia
using HiGHS
nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
mip_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver))
```

The feasibility pump is used at the start of Juniper to find a feasible solution before the branch and bound part starts.  For some classes of problems this can be a highly effective pre-processor.

# Citing Juniper

If you find Juniper useful in your work, we kindly request that you cite the following [paper](https://link.springer.com/chapter/10.1007/978-3-319-93031-2_27) or [technical report](https://arxiv.org/abs/1804.07332):

```
@inproceedings{juniper,
     Author = {Ole Kröger and Carleton Coffrin and Hassan Hijazi and Harsha Nagarajan},
     Title = {Juniper: An Open-Source Nonlinear Branch-and-Bound Solver in Julia},
     booktitle="Integration of Constraint Programming, Artificial Intelligence, and Operations Research",
     pages="377--386",
     year="2018",
     publisher="Springer International Publishing",
     isbn="978-3-319-93031-2"
}
```

# Juniper

[![CI](https://github.com/lanl-ansi/Juniper.jl/workflows/CI/badge.svg)](https://github.com/lanl-ansi/Juniper.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/lanl-ansi/Juniper.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/Juniper.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://lanl-ansi.github.io/Juniper.jl/stable)

[Juniper](http://github.com/lanl-ansi/Juniper.jl) (Jump Nonlinear Integer Program solver) is a solver for mixed-integer nonlinear programs. 

It is a heuristic which is not guaranteed to find the global optimum.
If you need the global optimum, check out [Alpine](http://github.com/lanl-ansi/Alpine.jl).

## Installation

Install Juniper using the Julia package manager:

```julia
import Pkg
Pkg.add("JuMP")
```

## Use with JuMP

Use Juniper with JuMP as follows:
```julia
using JuMP, Juniper, Ipopt
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt)
model = Model(optimizer)
v = [10, 20, 12, 23, 42]
w = [12, 45, 12, 22, 21]
@variable(model, x[1:5], Bin)
@objective(model, Max, v' * x)
@constraint(model, sum(w[i]*x[i]^2 for i in 1:5) <= 45)
optimize!(model)
println(termination_status(model))
println(objective_value(model))
println(value.(x))
```

The `nl_solver` is used by Juniper to solve continuous nonlinear sub-problems while Juniper searches for acceptable assignments to the discrete variables.
A common choice is [Ipopt](https://github.com/jump-dev/Ipopt.jl), but any optimizer that supports the continuous relaxation of the model may be used.

To solve problems with more complex nonlinear functions, use the `@NLconstraint` and `@NLobjective` JuMP macros.

## Feasibility pump

If Juniper has difficulty finding feasible solutions on your model, try adding a solver that supports integer variables (for example, [HiGHS](https://github.com/jump-dev/HiGHS.jl)) to run a _feasiblity pump_:

```julia
using JuMP, Juniper, Ipopt, HiGHS
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
highs = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
model = Model(
    optimizer_with_attributes(
        Juniper.Optimizer,
        "nl_solver" => ipopt,
        "mip_solver" => highs,
    ),
)
```

The feasibility pump is used at the start of Juniper to find a feasible solution before the branch and bound part starts. 
For some classes of problems this can be a highly effective pre-processor.

## Citing Juniper

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

# Juniper

Release: [![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://lanl-ansi.github.io/Juniper.jl/stable)

Dev: [![Build Status](https://travis-ci.org/lanl-ansi/Juniper.jl.svg?branch=master)](https://travis-ci.org/lanl-ansi/Juniper.jl) [![codecov](https://codecov.io/gh/lanl-ansi/Juniper.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/Juniper.jl)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://lanl-ansi.github.io/Juniper.jl/latest)

# Idea

You have a non linear problem with discrete variables (MINLP) and want some more control over the branch and bound part.
The relaxation should be solveable by any solver you prefer. Some solvers might not be able to solve the mixed integer part by themselves.

Juniper (Jump Nonlinear Integer Program solver) is a heuristic for non convex problems.
You need the global optimum? Check out [Alpine.jl](http://github.com/lanl-ansi/Alpine.jl)

# Basic usage

Juniper can be installed via:

`] add Juniper` for Julia v1.

Then adding it to your project by

`using Juniper`

You also have to import your NLP solver i.e.

`using Ipopt`

as well as [JuMP](http://www.juliaopt.org/JuMP.jl)

Define `Juniper` as the optimizer:

```
optimizer = Juniper.Optimizer
nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
```

And give it a go:

```
using LinearAlgebra # for the dot product
m = Model(optimizer_with_attributes(optimizer, "nl_solver"=>nl_solver))

v = [10,20,12,23,42]
w = [12,45,12,22,21]
@variable(m, x[1:5], Bin)

@objective(m, Max, dot(v,x))

@NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

optimize!(m)

# retrieve the objective value, corresponding x values and the status
println(JuMP.value.(x))
println(JuMP.objective_value(m))
println(JuMP.termination_status(m))
```

This solver is a NLP solver therefore you should have at least one `NLconstraint` or `NLobjective`.

It is recommended to specify a mip solver as well i.e.

```
using Cbc
optimizer = Juniper.Optimizer
nl_solver= optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
mip_solver = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
m = Model(optimizer_with_attributes(optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver))
```

Then the feasibility pump is used to find a feasible solution before the branch and bound part starts. This turned out to be highly effective.

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

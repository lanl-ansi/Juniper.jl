# Juniper

Dev: [![Build Status](https://travis-ci.org/lanl-ansi/Juniper.jl.svg?branch=master)](https://travis-ci.org/lanl-ansi/Juniper.jl) [![codecov](https://codecov.io/gh/lanl-ansi/Juniper.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/Juniper.jl)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://lanl-ansi.github.io/Juniper.jl/latest)


# Idea

You have a non linear problem with discrete variables (MINLP) and want some more control over the branch and bound part.
The relaxation should be solveable by any solver you prefer. Some solvers might not be able to solve the mixed integer part by themselves.

Juniper (Jump Nonlinear Integer Program solver) is a heuristic for non convex problems.
You need the global optimum? Check out [POD.jl](http://github.com/lanl-ansi/POD.jl)

# Basic usage

Version v0.1.0 can be installed via:

`Pkg.add("Juniper")`

Then adding it to your project by

`using Juniper`

You also have to import your NLP solver i.e.

`using Ipopt`

as well as [JuMP](http://www.juliaopt.org/JuMP.jl)

Define `JuniperSolver` as your solver:

```
solver = JuniperSolver(IpoptSolver(print_level=0))
```

And give it a go:

```
m = Model(solver=solver)

v = [10,20,12,23,42]
w = [12,45,12,22,21]
@variable(m, x[1:5], Bin)

@objective(m, Max, dot(v,x))

@NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

status = solve(m)
```

This solver is a NLP solver therefore you should have at least one `NLconstraint` or `NLobjective`.

It is recommended to specify a mip solver as well i.e.

```
using Cbc
solver = JuniperSolver(IpoptSolver(print_level=0);                       mip_solver=CbcSolver())
```

Then the feasibility pump is used to find a feasible solution before the branch and bound part starts. This turned out to be highly effective.
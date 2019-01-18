# Juniper

Juniper (Jump Non linear Integer Program solver) is a solver for MixedIntegerNonLinearPrograms (MINLPs) written in Julia.
Juniper solves these kind of problems using a NLP solver and then branch and bound. If the NLP solver isn't global optimal then Juniper is a heuristic. 
You need the global optimum? Check out [Alpine.jl](http://github.com/lanl-ansi/Alpine.jl)

# Why?
You have a non linear problem with discrete variables (MINLP) and want some more control over the branch and bound part? => Use Juniper

You have a really good solver for the relaxation and just want to solve problems with discrete variables as well? Just combine your solver with Juniper.

# Basic usage

The latest version can be installed via:

`Pkg.add("Juniper")`

or for Julia v0.7 and v1:

`] add Juniper` as `]` is used to get to interact with the package manager.

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

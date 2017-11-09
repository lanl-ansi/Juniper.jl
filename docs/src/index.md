# MINLPBnB

MINLPBnB is a solver for MixedIntegerNonLinearPrograms (MINLPs) written in Julia.
MINLPBnB solves these kind of problems using a NLP solver and then branch and bound.

# Why?
You have a non linear problem with discrete variables (MINLP) and want some more control over the branch and bound part? => Use MINLPBnB

You have a really good solver for the relaxation and just want to solve problems with discrete variables as well? Just combine your solver with MINLPBnB.

# Basic usage

It is currently not registered therefore you have to clone the package to be able to use it.

`Pkg.clone("http://github.com/Wikunia/MINLPBnB")`

Then adding it to your project by

`using MINLPBnB`

You also have to import your NLP solver i.e.

`using Ipopt`

as well as [JuMP](http://www.juliaopt.org/JuMP.jl)

Define `MINLPBnBSolver` as your solver:

```
solver = MINLPBnBSolver(IpoptSolver(print_level=0))
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

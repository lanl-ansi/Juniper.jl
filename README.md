# MINLPBnB

Dev: [![Build Status](https://travis-ci.org/Wikunia/MINLPBnB.svg?branch=master)](https://travis-ci.org/Wikunia/MINLPBnB) [![codecov](https://codecov.io/gh/Wikunia/MINLPBnB/branch/master/graph/badge.svg)](https://codecov.io/gh/Wikunia/MINLPBnB)

# Idea

You have a non linear problem with discrete variables (MINLP) and want some more control over the branch and bound part. 
The relaxation should be solveable by any solver you prefer. Some solvers might not be able to solve the mixed integer part by themselves.

# Basic usage

It is currently not registered therefore you have to clone the package to be able to use it.

`Pkg.clone("http://github.com/Wikunia/MINLPBnB")`

Then adding it to your project by

`using MINLPBnB`

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

As this solver is a NLP solver you should have at least one `NLconstraint` or `NLobjective`.


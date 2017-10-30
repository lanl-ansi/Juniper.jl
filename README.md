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

# Configuration

```
MINLPBnBSolver(IpoptSolver(print_level=0))
```

This is the most basic (non existent) configuration of the solver.

You can add options like doing the following:

```
minlpbnb = MINLPBnBSolver(IpoptSolver(print_level=0);
    branch_strategy=:StrongPseudoCost
)
```

In that example the strategy used for branching is defined.

In the following the options are explained. The type for the option is given after `::` and the default value in `[]`.

**Attention:**
The default values might change in the future after several tests were executed to determine the best overall options. 

## branch_strategy::Symbol [:StrongPseudoCost]

Possible values:

* `:MostInfeasible`
    * Branch on variables closest to 0.5
* `:PseudoCost`
    * Use `:MostInfeasible` first and then [Pseudo Cost Branching](https://en.wikipedia.org/wiki/Branch_and_cut#Branching_Strategies).
* `:StrongPseudoCost`
    * Use [Strong Branching](https://en.wikipedia.org/wiki/Branch_and_cut#Branching_Strategies) first and then `:PseudoCost`.

## Options for strong branching

### strong_branching_nvars::Int64 [5]

Defines the number of variables to consider for strong branching. 

### strong_branching_nsteps::Int64 [1]

Defines the number of steps in which strong branching is used. `:PseudoCost` will be used for later steps.

### strong_restart::Bool [true]

If a child while running strong branching is infeasible this holds for the whole node. Therefore we can tighten the bounds and rerun the strong branch part. (This might occur more then once)


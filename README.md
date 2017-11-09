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

# Configuration

```
MINLPBnBSolver(IpoptSolver(print_level=0))
```

This is the most basic configuration of the solver.

The first argument defines the solver for the relaxation here `IpoptSolver`. I use [Ipopt](https://projects.coin-or.org/Ipopt) for all the tests as well. The Ipopt julia package is described [here](https://github.com/JuliaOpt/Ipopt.jl) The solver itself can have parameters i.e `print_level=0`.

A list of some NLP solvers is mentioned [here](http://www.juliaopt.org/JuMP.jl/0.18/installation.html#getting-solvers)

You can add options doing the following:

```
minlpbnb = MINLPBnBSolver(IpoptSolver(print_level=0);
    branch_strategy=:StrongPseudoCost
)
```

In this example the strategy used for branching is defined.

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

## traverse_strategy::Symbol [:BFS]

Determines which node is used for the next branching.

Possible values:

* `:BFS`
    * Best-first-search: always take the node with the best bound
* `:DFS`
    * Depth-first-search: always take the node with the highest level
    * Might find a feasible solution faster
* `:DBFS`
    * Use of DFS first until the first feasible solution is found then switch to BFS

## incumbent_constr::Bool [true]

Add a constraint `objval >=/<= incumbent`. 

## Options for strong branching

### strong_branching_nvars::Int64 [5]

Defines the number of variables to consider for strong branching. 

### strong_branching_nsteps::Int64 [1]

Defines the number of steps in which strong branching is used. `:PseudoCost` will be used for later steps.

### strong_restart::Bool [true]

If a child while running strong branching is infeasible this holds for the whole node. Therefore we can tighten the bounds and rerun the strong branch part. (This might occur more then once)

## Parallel

MINLPBnB can be run in parallel to speed up the algorithm.
You have to start julia with `julia -p P` where `P` is the number of processors available or at least the number of processors you want to use.

Then you have to specify the number of processor as an option.

### processors::Int64 [1]

The number of processors used for the branch and bound part. **Attention:** Even if you start julia using
`julia -p P` you still have to define the number of processors using this option.

## Other options:

[UserLimits](readme/user_limits.md)

## Statistics 

[Statistics](readme/statistics.md)
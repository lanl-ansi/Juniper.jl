## General

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

## Branching

### branch_strategy::Symbol [:StrongPseudoCost]

Possible values:

* `:MostInfeasible`
    * Branch on variables closest to 0.5
* `:PseudoCost`
    * Use `:MostInfeasible` first and then [Pseudo Cost Branching](https://en.wikipedia.org/wiki/Branch_and_cut#Branching_Strategies).
* `:StrongPseudoCost`
    * Use [Strong Branching](https://en.wikipedia.org/wiki/Branch_and_cut#Branching_Strategies) first and then `:PseudoCost`.

### traverse_strategy::Symbol [:BFS]

Determines which node is used for the next branching.

Possible values:

* `:BFS`
    * Best-first-search: always take the node with the best bound
* `:DFS`
    * Depth-first-search: always take the node with the highest level
    * Might find a feasible solution faster
* `:DBFS`
    * Use of DFS first until the first feasible solution is found then switch to BFS

## Objective Cuts
### incumbent_constr::Bool [true]

Add a constraint `objval >=/<= incumbent`. 

### obj_epsilon::0 [Float64]

Add a constraint of the following form at the root node if not `0`.

If minimizing:

$\text{obj } \leq (1+\epsilon)\text{LB}$

If maximizing:

$\text{obj } \geq (1-\epsilon)\text{UB}$

## Options for strong branching

### strong_branching_perc::Float64 [25]

Defines the percentage of variables to consider for strong branching. 
If the number of variables is smaller than `2` it is fixed at `2` as strong branching doesn't make sense for one variable. 

### strong_branching_nsteps::Int64 [1]

Defines the number of steps in which strong branching is used. `:PseudoCost` will be used for later steps.

### strong_restart::Bool [true]

If a child, while running strong branching, is infeasible this holds for the whole node. Therefore we can tighten the bounds and rerun the strong branch part. (This might occur more then once)

## Options gain computation

The objective gain is used for pseudo cost and therefore also for strong pseudo cost branching.
$\mu$ is a parameter inside the computation of a so called score.
The version we use for the score function is described [here](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf).

### gain_mu::Float64 [0.167]

The parameter is used and a bit described in this [paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf) in (3).

## Parallel

MINLPBnB can be run in parallel to speed up the algorithm.
You have to start julia with `julia -p P` where `P` is the number of processors available or at least the number of processors you want to use.

Then you have to specify the number of processor as an option.

### processors::Int64 [1]

The number of processors used for the branch and bound part. **Attention:** Even if you start julia using
`julia -p P` you still have to define the number of processors using this option.

## User Limits

You can stop the solver before the optimal solution is found.
This is reasonable if the problem is to big to solve to optimality fast.
If the solver stops because of one of the following options the status `:UserLimit` is returned.

### time_limit::Float64 [Inf]

The maximum time in seconds the solver should run. 

**Note:** The solver will check after each branching step therefore it isn't a strict limit and depends on the duration of a relaxation.

### mip_gap::Float64 [0.01]%

Stops the solver if the gap is smaller than `mip_gap`. This is a percentage value.

### best_obj_stop::Float [NaN]

If an incumbent is found which is better than `best_obj_stop` the incumbent is returned. A warning gets thrown if `best_obj_stop` can't be reached.

### solution_limit::Int [0]

The solver stops if the requested amount of feasible solutions is found.
If `0` the option gets ignored.

## Logging

### log_levels::Vector{Symbol} [[:Table,:Info]]

You can change the option `log_levels` to define what kind of logs you want to see.

The default is `[:Table,:Info]` which outputs something like this:

![default-logging](https://user-images.githubusercontent.com/4931746/32625934-07b7db3c-c58e-11e7-922d-18a0a8776437.png)

Possible symbols which can be added to the vector are:

- :Timing
    - Provides some more timing informations
- :Options
    - prints all options that differ from default
- :AllOptions
    - prints all options 
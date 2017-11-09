## Options for user limits

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
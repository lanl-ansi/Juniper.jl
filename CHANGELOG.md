# Juniper.jl Changelog

### v0.9.2
- Bugfix: `DualStatus` now correctly returns `MOI.NO_SOLUTION` [PR #250](https://github.com/lanl-ansi/Juniper.jl/pull/250)
- Bugfix: fix setting default values for `MOI.TimeLimitSec` [PR #261](https://github.com/lanl-ansi/Juniper.jl/pull/261)
- Bugfix: FPump failed if `MOI.ScalarNonlinearFunction`s were present [PR #264](https://github.com/lanl-ansi/Juniper.jl/pull/264)

### v0.9.1
- Updated to JuMP v1.0

### v0.9.0
- Updated to MathOptInterface v1.0

### v0.8.0
- Updated to MathOptInterface v0.10

### v0.7.1
- Bugfix `set_subsolver_option!` nothing should be returned if `subsolver_name` doesn't match [PR #214](https://github.com/lanl-ansi/Juniper.jl/pull/214)

### v0.7.0
- Removed JuMP dependency [PR #207 by blegat](https://github.com/lanl-ansi/Juniper.jl/pull/207)

### v0.6.5
- Bugfix: FPump failed for an optimizer without attributes [Issue #201](https://github.com/lanl-ansi/Juniper.jl/pull/201)

### v0.6.4
- Use of `JUNIPER_RNG` for random numbers. Makes sure that `rand` outside of Juniper is unaffected.
- Option to set the random seed with the option `seed`.
- Bugfix: Generating the real nlp in the feasibility pump failed when there was no objective. [Issue #195](https://github.com/lanl-ansi/Juniper.jl/pull/195)

### 0.6.3
- Bugfix: `registered_functions` failed in feasibility pump [Issue #189](https://github.com/lanl-ansi/Juniper.jl/pull/189)

### 0.6.2
- Bugfix: `@NLexpression` failed in feasibility pump [Issue #184](https://github.com/lanl-ansi/Juniper.jl/pull/184)

### 0.6.1
- Bugfix: Use feasibility pump by default if MIP solver is used [#183](https://github.com/lanl-ansi/Juniper.jl/pull/183)

### 0.6.0
- Support for JuMP v0.21.0

### 0.5.3
- Documentation and useful error message for user defined functions
- Ability to set incumbent as `start_value`
- Support for NumberOfTreads & ResultCount
- Support for Interval constraints [#168](https://github.com/lanl-ansi/Juniper.jl/issues/168)

### 0.5.2
- bugfix: dereferencing when dereferencing twice (i.e. `@NLexpression`)

### 0.5.1
- Changed option `strong_branching_approx_time_limit` to `strong_branching_time_limit`
  - Default is still 100s

### 0.5.0
- Upgraded to MOI 0.9.1 and JuMP v0.20.0
  - Implementation of `TimeLimitSec`, `Silent` & `RawParameter`

### v0.4.3
- Bugfix: No bounds on binary variables in Feasibility pump [#143](https://github.com/lanl-ansi/Juniper.jl/issues/143)
- Bugfix: Feasibility pump if no objective in problem [#145](https://github.com/lanl-ansi/Juniper.jl/issues/145)

### v0.4.2
- Ability to not accept `almost` MOI solver status codes in tree search
- Better handling of MOI solver status codes internally
- Bugfix: Don't call objective value if sub problem not solved [#130](https://github.com/lanl-ansi/Juniper.jl/issues/130)
- Return ALMOST_LOCALLY_SOLVED if corresponding relaxation is only almost solved

### v0.4.1
- Support for user defined functions see issue [#118](https://github.com/lanl-ansi/Juniper.jl/issues/118)

### v0.4.0
- MPB -> MOI
- JuMP 0.19
- Dropped support for Julia versions < 1.0

### v0.3.0
- removes support for Julia < v1

### v0.2.6
- bugfix in init_strong_restart
- bugfix in mip_gap if objval=0
- bugfix if gap was 0 in table printing
- support for primal start values

### v0.2.5
- add linear constraints using @constraint in root model
- in BFS mode use the highest depth if several nodes have the same objective

### v0.2.4
- Remove Julia version upper bound

### v0.2.3
- Again support for `mu_init` in Ipopt

### v0.2.2
- Add support for Julia v1.0

### v0.2.1
- Add support for Julia v0.6 and v0.7

### v0.2.0
- Bugfix: Reset of `mu_init` in Ipopt options to have the default `mu_init` if `solve` is called again
- Bugfix: Break on time limit in relaxation, fpump and strong branching
- Bugfix: Infeasible in Reliability branching
- Strong branching:
    - Change bounds even if no restart
    - branch on best variable with two children

### v0.1.1
- Freemodel for commerical nlp solvers with license restrictions
- More convenient parallel options
    - `processors = 2` now uses 2 processors for solving nodes and one thread for supervision

### v0.1.0
- Traverse options (BFS,DFS,DBFS)
- Branch options (Strong, Reliable, Pseudo, MostInfeasible)
- Parallel solving of nodes
- Basic feasibility pump

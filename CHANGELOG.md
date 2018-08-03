# Juniper.jl Changelog

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

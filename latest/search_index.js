var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Juniper-1",
    "page": "Home",
    "title": "Juniper",
    "category": "section",
    "text": "Juniper (Jump Non linear Integer Program solver) is a solver for MixedIntegerNonLinearPrograms (MINLPs) written in Julia. Juniper solves these kind of problems using a NLP solver and then branch and bound. If the NLP solver isn't global optimal then Juniper is a heuristic.  You need the global optimum? Check out POD.jl"
},

{
    "location": "index.html#Why?-1",
    "page": "Home",
    "title": "Why?",
    "category": "section",
    "text": "You have a non linear problem with discrete variables (MINLP) and want some more control over the branch and bound part? => Use JuniperYou have a really good solver for the relaxation and just want to solve problems with discrete variables as well? Just combine your solver with Juniper."
},

{
    "location": "index.html#Basic-usage-1",
    "page": "Home",
    "title": "Basic usage",
    "category": "section",
    "text": "It is currently not registered therefore you have to clone the package to be able to use it.Pkg.clone(\"http://github.com/lanl-ansi/Juniper\")Then adding it to your project byusing JuniperYou also have to import your NLP solver i.e.using Ipoptas well as JuMPDefine JuniperSolver as your solver:solver = JuniperSolver(IpoptSolver(print_level=0))And give it a go:m = Model(solver=solver)\n\nv = [10,20,12,23,42]\nw = [12,45,12,22,21]\n@variable(m, x[1:5], Bin)\n\n@objective(m, Max, dot(v,x))\n\n@NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   \n\nstatus = solve(m)\nThis solver is a NLP solver therefore you should have at least one NLconstraint or NLobjective."
},

{
    "location": "options.html#",
    "page": "Options",
    "title": "Options",
    "category": "page",
    "text": ""
},

{
    "location": "options.html#General-1",
    "page": "Options",
    "title": "General",
    "category": "section",
    "text": "JuniperSolver(IpoptSolver(print_level=0))This is the most basic configuration of the solver.The first argument defines the solver for the relaxation here IpoptSolver. I use Ipopt for all the tests as well. The Ipopt julia package is described here The solver itself can have parameters i.e print_level=0.A list of some NLP solvers is mentioned hereYou can add options doing the following:juniper = JuniperSolver(IpoptSolver(print_level=0);\n    branch_strategy=:StrongPseudoCost\n)In this example the strategy used for branching is defined.In the following the options are explained. The type for the option is given after :: and the default value in [].Attention: The default values might change in the future after several tests were executed to determine the best overall options. "
},

{
    "location": "options.html#Branching-1",
    "page": "Options",
    "title": "Branching",
    "category": "section",
    "text": ""
},

{
    "location": "options.html#branch_strategy::Symbol-[:StrongPseudoCost]-1",
    "page": "Options",
    "title": "branch_strategy::Symbol [:StrongPseudoCost]",
    "category": "section",
    "text": "Possible values::MostInfeasible\nBranch on variables closest to 0.5\n:PseudoCost\nUse :MostInfeasible first and then Pseudo Cost Branching.\n:StrongPseudoCost\nUse Strong Branching first and then :PseudoCost.\n:Reliability\nUse Reliability Branching in a slightly different version."
},

{
    "location": "options.html#traverse_strategy::Symbol-[:BFS]-1",
    "page": "Options",
    "title": "traverse_strategy::Symbol [:BFS]",
    "category": "section",
    "text": "Determines which node is used for the next branching.Possible values::BFS\nBest-first-search: always take the node with the best bound\n:DFS\nDepth-first-search: always take the node with the highest level\nMight find a feasible solution faster\n:DBFS\nUse of DFS first until the first feasible solution is found then switch to BFS"
},

{
    "location": "options.html#Objective-Cuts-1",
    "page": "Options",
    "title": "Objective Cuts",
    "category": "section",
    "text": ""
},

{
    "location": "options.html#incumbent_constr::Bool-[true]-1",
    "page": "Options",
    "title": "incumbent_constr::Bool [true]",
    "category": "section",
    "text": "Add a constraint objval >=/<= incumbent. "
},

{
    "location": "options.html#obj_epsilon::0-[Float64]-1",
    "page": "Options",
    "title": "obj_epsilon::0 [Float64]",
    "category": "section",
    "text": "Add a constraint of the following form at the root node if not 0.If minimizing:textobj  leq (1+epsilon)textLBIf maximizing:textobj  geq (1-epsilon)textUB"
},

{
    "location": "options.html#Options-for-strong-branching-1",
    "page": "Options",
    "title": "Options for strong branching",
    "category": "section",
    "text": ""
},

{
    "location": "options.html#strong_branching_perc::Float64-[25]-1",
    "page": "Options",
    "title": "strong_branching_perc::Float64 [25]",
    "category": "section",
    "text": "Defines the percentage of variables to consider for strong branching.  If set to 25 it means that strong branching is performed on 25% of all discrete variables. Variables which are discrete in the relaxation aren't considered again but count to the number of  all discrete variables. If the number of variables is smaller than 2 it is fixed at 2 as strong branching doesn't make sense for one variable. "
},

{
    "location": "options.html#strong_branching_nsteps::Int64-[1]-1",
    "page": "Options",
    "title": "strong_branching_nsteps::Int64 [1]",
    "category": "section",
    "text": "Defines the number of steps in which strong branching is used. :PseudoCost will be used for later steps."
},

{
    "location": "options.html#strong_restart::Bool-[true]-1",
    "page": "Options",
    "title": "strong_restart::Bool [true]",
    "category": "section",
    "text": "If a child, while running strong branching, is infeasible this holds for the whole node. Therefore we can tighten the bounds and rerun the strong branch part. (This might occur more then once) This option is also used in reliablity branching."
},

{
    "location": "options.html#Options-for-reliablity-branching-1",
    "page": "Options",
    "title": "Options for reliablity branching",
    "category": "section",
    "text": "The implemented version of reliablity branching uses the gain score as in pseudo cost branching  and if some branching variables aren't reliable in a sense that strong branching wasn't performed  at least reliablility_branching_threshold times then strong branching is performed on those. Afterwards it will be branched on the variable with the highest gain score."
},

{
    "location": "options.html#reliablility_branching_perc::Float64-[25]-1",
    "page": "Options",
    "title": "reliablility_branching_perc::Float64 [25]",
    "category": "section",
    "text": "Defines the percentage of variables to consider for the strong branching part of reliablity branching. If the number of variables is smaller than 2 it is fixed at 2 as strong branching doesn't make sense for one variable. "
},

{
    "location": "options.html#reliablility_branching_threshold::Int64-[5]-1",
    "page": "Options",
    "title": "reliablility_branching_threshold::Int64 [5]",
    "category": "section",
    "text": "Defines whether strong branching is used to improve the reliability of the gain score. If a variable was used less than reliablility_branching_threshold times for strong branching then strong branching is performed on some of those candidates. The amount of candidates used is calculated by reliablility_branching_perc."
},

{
    "location": "options.html#strong_restart::Bool-[true]-2",
    "page": "Options",
    "title": "strong_restart::Bool [true]",
    "category": "section",
    "text": "This option is explained in strong branching but is also used for reliability branching."
},

{
    "location": "options.html#Options-gain-computation-1",
    "page": "Options",
    "title": "Options gain computation",
    "category": "section",
    "text": "The objective gain is used for pseudo cost and therefore also for strong pseudo cost branching. mu is a parameter inside the computation of a so called score. The version we use for the score function is described here."
},

{
    "location": "options.html#gain_mu::Float64-[0.167]-1",
    "page": "Options",
    "title": "gain_mu::Float64 [0.167]",
    "category": "section",
    "text": "The parameter is used and a bit described in this paper in (3)."
},

{
    "location": "options.html#Parallel-1",
    "page": "Options",
    "title": "Parallel",
    "category": "section",
    "text": "Juniper can be run in parallel to speed up the algorithm. You have to start julia with julia -p P where P is the number of processors available or at least the number of processors you want to use.Then you have to specify the number of processor as an option."
},

{
    "location": "options.html#processors::Int64-[1]-1",
    "page": "Options",
    "title": "processors::Int64 [1]",
    "category": "section",
    "text": "The number of processors used for the branch and bound part. Attention: Even if you start julia using julia -p P you still have to define the number of processors using this option."
},

{
    "location": "options.html#User-Limits-1",
    "page": "Options",
    "title": "User Limits",
    "category": "section",
    "text": "You can stop the solver before the optimal solution is found. This is reasonable if the problem is to big to solve to optimality fast. If the solver stops because of one of the following options the status :UserLimit is returned."
},

{
    "location": "options.html#time_limit::Float64-[Inf]-1",
    "page": "Options",
    "title": "time_limit::Float64 [Inf]",
    "category": "section",
    "text": "The maximum time in seconds the solver should run. Note: The solver will check after each branching step therefore it isn't a strict limit and depends on the duration of a relaxation."
},

{
    "location": "options.html#mip_gap::Float64-[0.0001]-1",
    "page": "Options",
    "title": "mip_gap::Float64 [0.0001]",
    "category": "section",
    "text": "Stops the solver if the gap is smaller than mip_gap. The default is 0.01%."
},

{
    "location": "options.html#best_obj_stop::Float-[NaN]-1",
    "page": "Options",
    "title": "best_obj_stop::Float [NaN]",
    "category": "section",
    "text": "If an incumbent is found which is better than best_obj_stop the incumbent is returned. A warning gets thrown if best_obj_stop can't be reached."
},

{
    "location": "options.html#solution_limit::Int-[0]-1",
    "page": "Options",
    "title": "solution_limit::Int [0]",
    "category": "section",
    "text": "The solver stops if the requested amount of feasible solutions is found. If 0 the option gets ignored."
},

{
    "location": "options.html#Logging-1",
    "page": "Options",
    "title": "Logging",
    "category": "section",
    "text": ""
},

{
    "location": "options.html#log_levels::Vector{Symbol}-[[:Table,:Info,:Options]]-1",
    "page": "Options",
    "title": "log_levels::Vector{Symbol} [[:Table,:Info,:Options]]",
    "category": "section",
    "text": "You can change the option log_levels to define what kind of logs you want to see.The output for [:Table,:Info] looks something like this:(Image: default-logging):Optionsincludes something like this before the info is printed:time_limit               : 10.0\nstrong_branching_nsteps  : 5Possible symbols which can be added to the vector are::Timing\nProvides some more timing informations\n:AllOptions\nprints all options "
},

{
    "location": "extras.html#",
    "page": "Extras",
    "title": "Extras",
    "category": "page",
    "text": ""
},

{
    "location": "extras.html#Statistics-1",
    "page": "Extras",
    "title": "Statistics",
    "category": "section",
    "text": "You can get some statistics after the problem is solved.Variable Description\nnintvars Number of integer variables\nnbinvars Number of binary variables\nnnodes Number of explored nodes in branch and bound\nncuts Number of cuts\nnbranches Number of branches\nnlevels Deepest level reached (Root node is level 1)To access these statistics you can use m.internalModel i.e.:m.internalModel.nintvars"
},

]}

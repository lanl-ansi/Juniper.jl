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
    "text": "Version v0.1.0 can be installed via:Pkg.add(\"Juniper\")Then adding it to your project byusing JuniperYou also have to import your NLP solver i.e.using Ipoptas well as JuMPDefine JuniperSolver as your solver:solver = JuniperSolver(IpoptSolver(print_level=0))And give it a go:m = Model(solver=solver)\n\nv = [10,20,12,23,42]\nw = [12,45,12,22,21]\n@variable(m, x[1:5], Bin)\n\n@objective(m, Max, dot(v,x))\n\n@NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   \n\nstatus = solve(m)\nThis solver is a NLP solver therefore you should have at least one NLconstraint or NLobjective."
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
    "text": "The most basic configuration of Juniper is:JuniperSolver(IpoptSolver(print_level=0))The first argument defines the solver for the relaxation here IpoptSolver. Ipopt is used for all our test cases. The Ipopt julia package is described here. The solver itself can have parameters i.e print_level=0.JuMP supports a lot of different NLP solvers (open source as well as commercial). A list of some NLP solvers is mentioned hereYou can add options doing the following:m = Model()\njuniper = JuniperSolver(IpoptSolver(print_level=0);\n    branch_strategy=:StrongPseudoCost\n)\nsetsolver(m, juniper)In this example the strategy used for branching is defined.In the following the options are explained. The type for the option is given after :: and the default value in [].Attention: The default values might change in the future after several tests were executed to determine the best overall options. "
},

{
    "location": "options.html#Tolerance-1",
    "page": "Options",
    "title": "Tolerance",
    "category": "section",
    "text": ""
},

{
    "location": "options.html#atol::Float64-[1e-6]-1",
    "page": "Options",
    "title": "atol::Float64 [1e-6]",
    "category": "section",
    "text": "This tolerance is used to check whether a value is integer or not and is used in the feasibility pump. More information about the feasibility pump can be found here."
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
    "text": "Possible values::MostInfeasible\nBranch on variables closest to 0.5\n:PseudoCost\nUse :MostInfeasible first and then Pseudo Cost Branching.\n:StrongPseudoCost\nUse Strong Branching first and then :PseudoCost\nMore options for strong branching are described here\n:Reliability\nUse Reliability Branching in a slightly different version."
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
    "location": "options.html#incumbent_constr::Bool-[false]-1",
    "page": "Options",
    "title": "incumbent_constr::Bool [false]",
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
    "location": "options.html#strong_branching_perc::Float64-[100]-1",
    "page": "Options",
    "title": "strong_branching_perc::Float64 [100]",
    "category": "section",
    "text": "Defines the percentage of variables to consider for strong branching.  If set to 25 it means that strong branching is performed on 25% of all discrete variables. Variables which are discrete in the relaxation aren't considered again but count to the number of  all discrete variables. If the number of variables is smaller than 2 it is fixed at 2 as strong branching doesn't make sense for one variable. Attention: strong_branching_approx_time_limit might change this value."
},

{
    "location": "options.html#strong_branching_nsteps::Int64-[1]-1",
    "page": "Options",
    "title": "strong_branching_nsteps::Int64 [1]",
    "category": "section",
    "text": "Defines the number of steps in which strong branching is used. :PseudoCost will be used for later steps."
},

{
    "location": "options.html#strong_branching_approx_time_limit::Float64-[100]s-1",
    "page": "Options",
    "title": "strong_branching_approx_time_limit::Float64 [100]s",
    "category": "section",
    "text": "For big problems with either a lot of variables or a long relaxation time it turned out to be reasonable to reduce the number of strong branching variables.A small example:strong_branching_perc is set to 100%.\nThe root relaxation takes 5 seconds and there are 100 discrete variables.\nNow the approximated time for strong branching (without considering restarts) is   2*5s*100 = 1000s because each discrete variable has two children.By using strong_branching_approx_time_limit = 100 the number of strong branching variables is reduced to 10 because 2*5s*10 = 100. If you don't want to use this time limit you can set it to Inf."
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
    "location": "options.html#Feasibility-Pump-1",
    "page": "Options",
    "title": "Feasibility Pump",
    "category": "section",
    "text": "Juniper has the option to find a feasible solution before the branch and bound part starts. The following options describe how to use the feasibility pump."
},

{
    "location": "options.html#feasibility_pump::Bool-[Auto]-1",
    "page": "Options",
    "title": "feasibility_pump::Bool [Auto]",
    "category": "section",
    "text": "Determines whether or not the feasibility pump should be used to get a feasible solution. The default is true if a mip solver is specified and false otherwise. Attention: If set to true you need to also set the mip_solver option."
},

{
    "location": "options.html#mip_solver::MathProgBase.AbstractMathProgSolver-[nothing]-1",
    "page": "Options",
    "title": "mip_solver::MathProgBase.AbstractMathProgSolver [nothing]",
    "category": "section",
    "text": "This has to be set to a mip solver if the feasibility pump should be used. A list of some MIP solvers is mentioned here.If you want to use GLPK you would need to useusing GLPKMathProgInterfaceand set the option with mip_solver=GLPKSolverMIP()"
},

{
    "location": "options.html#feasibility_pump_time_limit::Int64-[60]s-1",
    "page": "Options",
    "title": "feasibility_pump_time_limit::Int64 [60]s",
    "category": "section",
    "text": "The time limit of the feasibility pump in seconds. After that time limit the branch and bound part starts whether a feasible solution was found or not."
},

{
    "location": "options.html#feasibility_pump_tolerance_counter::Int64-[5]-1",
    "page": "Options",
    "title": "feasibility_pump_tolerance_counter::Int64 [5]",
    "category": "section",
    "text": "In the feasibility pump the objective is to reduce the difference between the mip and the nlp solution. If the default tolerance (atol) can't be reached for feasibility_pump_tolerance_counter consecutive times but the tolerance*10 can be reached. The tolerance will be switched after feasibility_pump_tolerance_counter and a warning will be thrown. If there is no warning like Real objective wasn't solved to optimality afterwards there is no need to worry at all. Normally in the end of the feasibility pump the real objective is used to improve the  objective value. If this is possible the warning before can be ignored as it is given that the solution has no rounding issues.  If this can't be done a warning like Real objective wasn't solved to optimality is thrown. This means that the objective might be not the best possible given the mip solution and if a warning for the tolerance change was thrown there might be rounding issues.  You can set this value to a huge number (more than 100 should be enough) if you don't want to use this option."
},

{
    "location": "options.html#tabu_list_length::Int64-[30]-1",
    "page": "Options",
    "title": "tabu_list_length::Int64 [30]",
    "category": "section",
    "text": "During the run of the feasibility pump it might happen that the alternating solve steps get into a cycle. By using a tabu list cycles can be avoided. The length determines the length of the cycle which will be avoided. If a cycle is encountered which is longer the feasibility pump terminates."
},

{
    "location": "options.html#num_resolve_nlp_feasibility_pump::Int64-[1]-1",
    "page": "Options",
    "title": "num_resolve_nlp_feasibility_pump::Int64 [1]",
    "category": "section",
    "text": "If the NLP is infeasible during the feasibility pump it can be restarted with a random starting point for the NL solver. This will be done as long as it is infeasible or num_resolve_nlp_feasibility_pump is reached."
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
    "location": "options.html#Resolve-1",
    "page": "Options",
    "title": "Resolve",
    "category": "section",
    "text": "Sometimes the non linear solver doesn't find a feasible solution in the first run."
},

{
    "location": "options.html#num_resolve_root_relaxation::Int-[3]-1",
    "page": "Options",
    "title": "num_resolve_root_relaxation::Int [3]",
    "category": "section",
    "text": "This especially bad if this happens for the root relaxation. If there is no optimal/local optimal solution in the root relaxation you can use this option to resolve a couple of time until a solution is found or the number of resolves exceeded this value."
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
    "text": "You can change the option log_levels to define what kind of logs you want to see.The output for [:Table,:Info] looks something like this:(Image: default-logging):Table#ONodes\nThe number of open nodes\nCLevel\nThe current node is at level ... of the tree\nIncumbent\nBest integral solution found\nBest Bound\nThe best bound of the open nodes\nGap \nThe gap between Incumbent and Best Bound\nTime\nThe time spend since the beginning of branch and bound\nDoesn't count time before branch and bound starts (i.e. feasibility pump or root relaxation)\nGainGap\nThe difference in percentage between a guessed gain and the actual gain.\nUsed if branch_strategy = PseudoCost or after strong branching / reliability branching.:Optionsincludes something like this before the info is printed:time_limit               : 10.0\nstrong_branching_nsteps  : 5Possible symbols which can be added to the vector are::Timing\nProvides some more timing information\n:AllOptions\nprints all options "
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

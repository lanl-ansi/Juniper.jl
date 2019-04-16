using JuMP, Juniper, Ipopt, Cbc
m = Model()

mutable struct blink
    cost::Float64
    x::Float64
    t::Float64
end

cost = randn(5)
x = randn(5)
t = randn(5)

Blinks = [blink(cost[i], x[i], t[i]) for i = 1:5]
p = [(1,3)]
@variables m begin
    δ[a in Blinks], Bin,(start = (a.cost == 0))
end

#Some constraints
for i in p
    @constraint(m, sum(δ[Blinks[j]] for j in i) <= 1)
end

#Budget constraint
B = 90
@constraint(m,Bud, sum(δ[a]*a.cost for a in Blinks) <= B)

function low_level(Blinks,δ...)
    ind = Int[]
    for i = 1:length(δ)
        if isapprox(δ[i], round(δ[i]))
            push!(ind,i)
        else
            Blinks[i].x = 0.0
            Blinks[i].t = 0.0
        end
    end

    #SolveEquilibrium(Nodes,vcat(Links,Clinks[ind]),ODs)
    #SolveEquilibrium modifies x and t for each member. That problem solves fine
    return  sum(a.x*a.t for a in Blinks)
end

#Function to be registered
ll(δ...) = low_level(Blinks,δ...)

JuMP.register(m,:ll,length(Blinks),ll,autodiff=true)
solver = JuniperSolver(IpoptSolver(print_level=0);
                       mip_solver=CbcSolver(),
                       registered_functions=[Juniper.register(:ll,length(Blinks),ll,autodiff=true)])

setsolver(m, solver)                

JuMP.setNLobjective(m,:Min,:($(Expr(:call, :ll, (δ[a] for a in Blinks)...))))

status = solve(m)
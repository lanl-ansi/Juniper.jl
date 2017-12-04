function generate_random_restart(m)
    values = []
    for i=1:m.num_var
        lbi_def = true
        ubi_def = true
        if m.l_var[i] > typemin(Int64) 
            lbi = m.l_var[i]
        else
            lbi = typemin(Int64)
            lbi_def = false
        end

        if m.l_var[i] < typemax(Int64) 
            ubi = m.l_var[i]
        else
            ubi = typemin(Int64)
            ubi_def = false
        end

        if !ubi_def && !lbi_def 
            ubi = 10
            lbi = -10
        elseif !ubi_def
            ubi = lbi+20
        elseif !lbi_def
            lbi = ubi-20
        end             

        if m.var_type[i] == :Cont
            push!(values,(ubi-lbi)*rand()+lbi)
        else
            push!(values, rand(lbi:ubi))
        end
    end
    return values
end
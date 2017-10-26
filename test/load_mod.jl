function model_from_mod(filename)
    f = open("./test/"*filename);
    for ln in eachline(f)
        ln = lstrip(ln)
        ln = rstrip(ln)
        if length(ln) == 0
            continue
        end
        if ln[1] == '#'
            continue
        end
        if ln[1:3] == "var"
            ln = replace(ln, r" :=\s*(\d*(\.?\d*)),","")
            ln = replace(ln, "var ","")
            ln = replace(ln, r"^(\S*)\s*>=\s*(\S*),\s*<=\s*(\S*);$",s"\2 <= \1 <= \3")      
            ln = replace(ln, r"^(\S*) binary.*",s"\1, Bin")            
            ln = "@variable(m, "*ln*")"
            println(ln)
        elseif ln[1] == 'e'
            ln = replace(ln, " = ", " == ")      
            ln = replace(ln, r"e\S*\s*","")  
            ln = ln[1:end-1] 
            ln = "@NLconstraint(m, "*ln*")"   
            println(ln)
        end
    end
end
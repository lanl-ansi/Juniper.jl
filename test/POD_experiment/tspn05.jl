function get_tspn05()

m = Model()

# ----- Variables ----- #
@variable(m, objvar)
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
@variable(m, x[x_Idx])
b_Idx = Any[11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
@variable(m, b[b_Idx], Bin)
JuMP.set_lower_bound(x[1], 68.0)
JuMP.set_upper_bound(x[1], 71.0)
JuMP.set_lower_bound(x[2], 65.0)
JuMP.set_upper_bound(x[2], 87.0)
JuMP.set_lower_bound(x[3], 107.0)
JuMP.set_upper_bound(x[3], 126.0)
JuMP.set_lower_bound(x[4], 38.0)
JuMP.set_upper_bound(x[4], 49.0)
JuMP.set_lower_bound(x[5], 40.0)
JuMP.set_upper_bound(x[5], 55.0)
JuMP.set_lower_bound(x[6], 54.0)
JuMP.set_upper_bound(x[6], 68.0)
JuMP.set_lower_bound(x[7], 92.0)
JuMP.set_upper_bound(x[7], 106.0)
JuMP.set_lower_bound(x[8], 113.0)
JuMP.set_upper_bound(x[8], 117.0)
JuMP.set_lower_bound(x[9], 82.0)
JuMP.set_upper_bound(x[9], 87.0)
JuMP.set_lower_bound(x[10], 76.0)
JuMP.set_upper_bound(x[10], 85.0)


# ----- Constraints ----- #
@NLconstraint(m, e1, sqrt( (x[1]-x[3])^2+ (x[2]-x[4])^2)*b[11]+sqrt( (x[1]-x[5])^2+ (x[2]-x[6])^2)*b[12]+sqrt( (x[1]-x[7])^2+ (x[2]-x[8])^2)*b[13]+sqrt( (x[1]-x[9])^2+ (x[2]-x[10])^2)*b[14]+sqrt( (x[3]-x[5])^2+ (x[4]-x[6])^2)*b[15]+sqrt( (x[3]-x[7])^2+ (x[4]-x[8])^2)*b[16]+sqrt( (x[3]-x[9])^2+ (x[4]-x[10])^2)*b[17]+sqrt( (x[5]-x[7])^2+ (x[6]-x[8])^2)*b[18]+sqrt( (x[5]-x[9])^2+ (x[6]-x[10])^2)*b[19]+sqrt( (x[7]-x[9])^2+ (x[8]-x[10])^2)*b[20]-objvar == 0.0)
@NLconstraint(m, e2, 0.444444444444444* (x[1])^2-61.7777777777778*x[1]+0.00826446280991736* (x[2])^2-1.25619834710744*x[2] <= -2193.51331496786)
@NLconstraint(m, e3, 0.0110803324099723* (x[3])^2-2.58171745152355*x[3]+0.0330578512396694* (x[4])^2-2.87603305785124*x[4] <= -211.938760559511)
@NLconstraint(m, e4, 0.0177777777777778* (x[5])^2-1.68888888888889*x[5]+0.0204081632653061* (x[6])^2-2.48979591836735*x[6] <= -115.049886621315)
@NLconstraint(m, e5, 0.0204081632653061* (x[7])^2-4.04081632653061*x[7]+0.25* (x[8])^2-57.5*x[8] <= -3505.27040816327)
@NLconstraint(m, e6, 0.16* (x[9])^2-27.04*x[9]+0.0493827160493827* (x[10])^2-7.95061728395062*x[10] <= -1461.45234567901)
@constraint(m, e7, b[11]+b[12]+b[13]+b[14] == 2.0)
@constraint(m, e8, b[11]+b[15]+b[16]+b[17] == 2.0)
@constraint(m, e9, b[12]+b[15]+b[18]+b[19] == 2.0)
@constraint(m, e10, b[13]+b[16]+b[18]+b[20] == 2.0)
@constraint(m, e11, b[14]+b[17]+b[19]+b[20] == 2.0)


# ----- Objective ----- #
@objective(m, Min, objvar)

return m
end
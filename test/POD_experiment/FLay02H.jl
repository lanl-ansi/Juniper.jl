function get_FLay02H()
    m = Model()

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
    ]
    @variable(m, x[x_Idx])
    b_Idx = Any[43, 44, 45, 46]
    @variable(m, b[b_Idx], Bin)
    JuMP.set_lower_bound(x[36], 0.0)
    JuMP.set_lower_bound(x[4], 0.0)
    JuMP.set_lower_bound(x[16], 0.0)
    JuMP.set_lower_bound(x[32], 0.0)
    JuMP.set_lower_bound(x[27], 0.0)
    JuMP.set_lower_bound(x[14], 0.0)
    JuMP.set_lower_bound(x[17], 0.0)
    JuMP.set_lower_bound(x[3], 0.0)
    JuMP.set_lower_bound(x[25], 0.0)
    JuMP.set_lower_bound(x[38], 0.0)
    JuMP.set_lower_bound(x[30], 0.0)
    JuMP.set_lower_bound(x[26], 0.0)
    JuMP.set_lower_bound(x[23], 0.0)
    JuMP.set_lower_bound(x[42], 0.0)
    JuMP.set_lower_bound(x[34], 0.0)
    JuMP.set_lower_bound(x[11], 0.0)
    JuMP.set_lower_bound(x[29], 0.0)
    JuMP.set_lower_bound(x[22], 0.0)
    JuMP.set_lower_bound(x[12], 0.0)
    JuMP.set_lower_bound(x[37], 0.0)
    JuMP.set_lower_bound(x[19], 0.0)
    JuMP.set_lower_bound(x[40], 0.0)
    JuMP.set_lower_bound(x[2], 0.0)
    JuMP.set_lower_bound(x[20], 0.0)
    JuMP.set_lower_bound(x[24], 0.0)
    JuMP.set_lower_bound(x[41], 0.0)
    JuMP.set_lower_bound(x[39], 0.0)
    JuMP.set_lower_bound(x[31], 0.0)
    JuMP.set_lower_bound(x[18], 0.0)
    JuMP.set_lower_bound(x[9], 0.0)
    JuMP.set_lower_bound(x[15], 0.0)
    JuMP.set_lower_bound(x[1], 0.0)
    JuMP.set_lower_bound(x[33], 0.0)
    JuMP.set_lower_bound(x[13], 0.0)
    JuMP.set_lower_bound(x[21], 0.0)
    JuMP.set_lower_bound(x[28], 0.0)
    JuMP.set_lower_bound(x[35], 0.0)
    JuMP.set_lower_bound(x[10], 0.0)
    JuMP.set_upper_bound(x[1], 29.0)
    JuMP.set_upper_bound(x[2], 29.0)
    JuMP.set_upper_bound(x[3], 29.0)
    JuMP.set_upper_bound(x[4], 29.0)
    JuMP.set_lower_bound(x[5], 1.0)
    JuMP.set_upper_bound(x[5], 40.0)
    JuMP.set_lower_bound(x[6], 1.0)
    JuMP.set_upper_bound(x[6], 50.0)
    JuMP.set_lower_bound(x[7], 1.0)
    JuMP.set_upper_bound(x[7], 40.0)
    JuMP.set_lower_bound(x[8], 1.0)
    JuMP.set_upper_bound(x[8], 50.0)
    JuMP.set_upper_bound(x[9], 30.0)
    JuMP.set_upper_bound(x[10], 30.0)

    # ----- Constraints ----- #
    @constraint(m, e1, -2 * x[9] - 2 * x[10] + objvar == 0.0)
    @constraint(m, e2, -x[1] - x[5] + x[9] >= 0.0)
    @constraint(m, e3, -x[2] - x[6] + x[9] >= 0.0)
    @constraint(m, e4, -x[3] - x[7] + x[10] >= 0.0)
    @constraint(m, e5, -x[4] - x[8] + x[10] >= 0.0)
    @NLconstraint(m, e6, 40 / x[7] - x[5] <= 0.0)
    @NLconstraint(m, e7, 50 / x[8] - x[6] <= 0.0)
    @constraint(m, e8, x[1] - x[11] - x[12] - x[13] - x[14] == 0.0)
    @constraint(m, e9, x[2] - x[15] - x[16] - x[17] - x[18] == 0.0)
    @constraint(m, e10, x[3] - x[19] - x[20] - x[21] - x[22] == 0.0)
    @constraint(m, e11, x[4] - x[23] - x[24] - x[25] - x[26] == 0.0)
    @constraint(m, e12, x[5] - x[27] - x[28] - x[29] - x[30] == 0.0)
    @constraint(m, e13, x[6] - x[31] - x[32] - x[33] - x[34] == 0.0)
    @constraint(m, e14, x[7] - x[35] - x[36] - x[37] - x[38] == 0.0)
    @constraint(m, e15, x[8] - x[39] - x[40] - x[41] - x[42] == 0.0)
    @constraint(m, e16, x[11] - 29 * b[43] <= 0.0)
    @constraint(m, e17, x[12] - 29 * b[44] <= 0.0)
    @constraint(m, e18, x[13] - 29 * b[45] <= 0.0)
    @constraint(m, e19, x[14] - 29 * b[46] <= 0.0)
    @constraint(m, e20, x[15] - 29 * b[43] <= 0.0)
    @constraint(m, e21, x[16] - 29 * b[44] <= 0.0)
    @constraint(m, e22, x[17] - 29 * b[45] <= 0.0)
    @constraint(m, e23, x[18] - 29 * b[46] <= 0.0)
    @constraint(m, e24, x[19] - 29 * b[43] <= 0.0)
    @constraint(m, e25, x[20] - 29 * b[44] <= 0.0)
    @constraint(m, e26, x[21] - 29 * b[45] <= 0.0)
    @constraint(m, e27, x[22] - 29 * b[46] <= 0.0)
    @constraint(m, e28, x[23] - 29 * b[43] <= 0.0)
    @constraint(m, e29, x[24] - 29 * b[44] <= 0.0)
    @constraint(m, e30, x[25] - 29 * b[45] <= 0.0)
    @constraint(m, e31, x[26] - 29 * b[46] <= 0.0)
    @constraint(m, e32, x[27] - 40 * b[43] <= 0.0)
    @constraint(m, e33, x[28] - 40 * b[44] <= 0.0)
    @constraint(m, e34, x[29] - 40 * b[45] <= 0.0)
    @constraint(m, e35, x[30] - 40 * b[46] <= 0.0)
    @constraint(m, e36, x[31] - 40 * b[43] <= 0.0)
    @constraint(m, e37, x[32] - 40 * b[44] <= 0.0)
    @constraint(m, e38, x[33] - 40 * b[45] <= 0.0)
    @constraint(m, e39, x[34] - 40 * b[46] <= 0.0)
    @constraint(m, e40, x[35] - 40 * b[43] <= 0.0)
    @constraint(m, e41, x[36] - 40 * b[44] <= 0.0)
    @constraint(m, e42, x[37] - 40 * b[45] <= 0.0)
    @constraint(m, e43, x[38] - 40 * b[46] <= 0.0)
    @constraint(m, e44, x[39] - 40 * b[43] <= 0.0)
    @constraint(m, e45, x[40] - 40 * b[44] <= 0.0)
    @constraint(m, e46, x[41] - 40 * b[45] <= 0.0)
    @constraint(m, e47, x[42] - 40 * b[46] <= 0.0)
    @constraint(m, e48, x[11] - x[15] + x[27] <= 0.0)
    @constraint(m, e49, -x[12] + x[16] + x[32] <= 0.0)
    @constraint(m, e50, x[21] - x[25] + x[37] <= 0.0)
    @constraint(m, e51, -x[22] + x[26] + x[42] <= 0.0)
    @constraint(m, e52, b[43] + b[44] + b[45] + b[46] == 1.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

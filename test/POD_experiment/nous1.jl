function get_nous1()
    m = Model()

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[
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
        43,
        44,
        45,
        46,
        47,
        48,
        49,
    ]
    @variable(m, x[x_Idx])
    b_Idx = Any[50, 51]
    @variable(m, b[b_Idx], Bin)
    JuMP.set_lower_bound(x[16], 0.0)
    JuMP.set_lower_bound(x[14], 0.0)
    JuMP.set_lower_bound(x[38], 0.0)
    JuMP.set_lower_bound(x[42], 0.0)
    JuMP.set_lower_bound(x[22], 0.0)
    JuMP.set_lower_bound(x[2], 0.0)
    JuMP.set_lower_bound(x[9], 0.0)
    JuMP.set_lower_bound(x[8], 0.0)
    JuMP.set_lower_bound(x[43], 0.0)
    JuMP.set_lower_bound(x[36], 0.0)
    JuMP.set_lower_bound(x[4], 0.0)
    JuMP.set_lower_bound(x[32], 0.0)
    JuMP.set_lower_bound(x[27], 0.0)
    JuMP.set_lower_bound(x[3], 0.0)
    JuMP.set_lower_bound(x[25], 0.0)
    JuMP.set_lower_bound(x[30], 0.0)
    JuMP.set_lower_bound(x[11], 0.0)
    JuMP.set_lower_bound(x[29], 0.0)
    JuMP.set_lower_bound(x[5], 0.0)
    JuMP.set_lower_bound(x[37], 0.0)
    JuMP.set_lower_bound(x[24], 0.0)
    JuMP.set_lower_bound(x[41], 0.0)
    JuMP.set_lower_bound(x[18], 0.0)
    JuMP.set_lower_bound(x[7], 0.0)
    JuMP.set_lower_bound(x[13], 0.0)
    JuMP.set_lower_bound(x[21], 0.0)
    JuMP.set_lower_bound(x[10], 0.0)
    JuMP.set_lower_bound(x[26], 0.0)
    JuMP.set_lower_bound(x[45], 0.0)
    JuMP.set_lower_bound(x[12], 0.0)
    JuMP.set_lower_bound(x[40], 0.0)
    JuMP.set_lower_bound(x[44], 0.0)
    JuMP.set_lower_bound(x[31], 0.0)
    JuMP.set_lower_bound(x[33], 0.0)
    JuMP.set_lower_bound(x[28], 0.0)
    JuMP.set_lower_bound(x[35], 0.0)
    JuMP.set_lower_bound(x[6], 0.0)
    JuMP.set_lower_bound(x[17], 0.0)
    JuMP.set_lower_bound(x[23], 0.0)
    JuMP.set_lower_bound(x[34], 0.0)
    JuMP.set_lower_bound(x[19], 0.0)
    JuMP.set_lower_bound(x[20], 0.0)
    JuMP.set_lower_bound(x[39], 0.0)
    JuMP.set_lower_bound(x[15], 0.0)
    JuMP.set_upper_bound(x[2], 300.0)
    JuMP.set_upper_bound(x[3], 300.0)
    JuMP.set_upper_bound(x[4], 300.0)
    JuMP.set_upper_bound(x[5], 300.0)
    JuMP.set_upper_bound(x[6], 300.0)
    JuMP.set_upper_bound(x[7], 300.0)
    JuMP.set_upper_bound(x[8], 300.0)
    JuMP.set_upper_bound(x[9], 300.0)
    JuMP.set_upper_bound(x[10], 300.0)
    JuMP.set_upper_bound(x[11], 300.0)
    JuMP.set_upper_bound(x[12], 300.0)
    JuMP.set_upper_bound(x[13], 300.0)
    JuMP.set_upper_bound(x[14], 300.0)
    JuMP.set_upper_bound(x[15], 300.0)
    JuMP.set_upper_bound(x[16], 300.0)
    JuMP.set_upper_bound(x[17], 300.0)
    JuMP.set_upper_bound(x[18], 300.0)
    JuMP.set_upper_bound(x[19], 300.0)
    JuMP.set_upper_bound(x[22], 100.0)
    JuMP.set_upper_bound(x[23], 100.0)
    JuMP.set_upper_bound(x[24], 100.0)
    JuMP.set_upper_bound(x[25], 100.0)
    JuMP.set_upper_bound(x[26], 100.0)
    JuMP.set_upper_bound(x[27], 100.0)
    JuMP.set_upper_bound(x[28], 1.0)
    JuMP.set_upper_bound(x[29], 1.0)
    JuMP.set_upper_bound(x[30], 1.0)
    JuMP.set_upper_bound(x[31], 1.0)
    JuMP.set_upper_bound(x[32], 1.0)
    JuMP.set_upper_bound(x[33], 1.0)
    JuMP.set_upper_bound(x[34], 1.0)
    JuMP.set_upper_bound(x[35], 1.0)
    JuMP.set_upper_bound(x[36], 1.0)
    JuMP.set_upper_bound(x[37], 1.0)
    JuMP.set_upper_bound(x[38], 1.0)
    JuMP.set_upper_bound(x[39], 1.0)
    JuMP.set_upper_bound(x[40], 1.0)
    JuMP.set_upper_bound(x[41], 1.0)
    JuMP.set_upper_bound(x[42], 1.0)
    JuMP.set_upper_bound(x[43], 1.0)
    JuMP.set_upper_bound(x[44], 1.0)
    JuMP.set_upper_bound(x[45], 1.0)
    JuMP.set_lower_bound(x[46], 0.85)
    JuMP.set_upper_bound(x[46], 1.0)
    JuMP.set_lower_bound(x[47], 0.85)
    JuMP.set_upper_bound(x[47], 1.0)
    JuMP.set_lower_bound(x[48], 0.85)
    JuMP.set_upper_bound(x[48], 1.0)
    JuMP.set_lower_bound(x[49], 0.85)
    JuMP.set_upper_bound(x[49], 1.0)

    # ----- Constraints ----- #
    @NLconstraint(
        m,
        e1,
        -(
            (
                0.0042656 * x[29] - 0.0005719 * x[28] +
                0.0093514 * x[46] +
                0.0077308 * x[47] - 0.0139904
            ) * x[4] +
            (
                0.0016371 * x[31] +
                0.0288996 * x[32] +
                0.0338147 * x[48] +
                0.0373349 * x[49] - 0.0661588
            ) * x[5]
        ) + objvar - 0.23947 * b[50] - 0.75835 * b[51] == 0.0
    )
    @constraint(m, e2, x[2] + x[3] + x[20] + x[21] == 300.0)
    @constraint(m, e3, x[6] - x[12] - x[13] == 0.0)
    @constraint(m, e4, x[7] - x[11] - x[14] - x[15] == 0.0)
    @constraint(m, e5, x[8] - x[10] - x[16] - x[17] == 0.0)
    @constraint(m, e6, x[9] - x[18] - x[19] == 0.0)
    @NLconstraint(
        m,
        e7,
        -x[10] * x[40] - 0.333333333333333 * x[2] + x[22] == 0.0
    )
    @NLconstraint(
        m,
        e8,
        -x[10] * x[41] - 0.333333333333333 * x[2] + x[23] == 0.0
    )
    @NLconstraint(
        m,
        e9,
        -x[10] * x[42] - 0.333333333333333 * x[2] + x[24] == 0.0
    )
    @NLconstraint(
        m,
        e10,
        -x[11] * x[37] - 0.333333333333333 * x[3] + x[25] == 0.0
    )
    @NLconstraint(
        m,
        e11,
        -x[11] * x[38] - 0.333333333333333 * x[3] + x[26] == 0.0
    )
    @NLconstraint(
        m,
        e12,
        -x[11] * x[39] - 0.333333333333333 * x[3] + x[27] == 0.0
    )
    @NLconstraint(m, e13, -(x[6] * x[34] + x[7] * x[37]) + x[22] == 0.0)
    @NLconstraint(m, e14, -(x[6] * x[35] + x[7] * x[38]) + x[23] == 0.0)
    @NLconstraint(m, e15, -(x[6] * x[36] + x[7] * x[39]) + x[24] == 0.0)
    @NLconstraint(m, e16, -(x[8] * x[40] + x[9] * x[43]) + x[25] == 0.0)
    @NLconstraint(m, e17, -(x[8] * x[41] + x[9] * x[44]) + x[26] == 0.0)
    @NLconstraint(m, e18, -(x[8] * x[42] + x[9] * x[45]) + x[27] == 0.0)
    @NLconstraint(m, e19, x[22] * x[46] - x[6] * x[34] == 0.0)
    @NLconstraint(m, e20, x[23] * x[47] - x[7] * x[38] == 0.0)
    @NLconstraint(m, e21, x[26] * x[48] - x[8] * x[41] == 0.0)
    @NLconstraint(m, e22, x[27] * x[49] - x[9] * x[45] == 0.0)
    @NLconstraint(
        m,
        e23,
        x[12] * x[34] +
        x[14] * x[37] +
        x[16] * x[40] +
        x[18] * x[43] +
        0.333333333333333 * x[20] == 30.0
    )
    @NLconstraint(
        m,
        e24,
        x[12] * x[35] +
        x[14] * x[38] +
        x[16] * x[41] +
        x[18] * x[44] +
        0.333333333333333 * x[20] == 50.0
    )
    @NLconstraint(
        m,
        e25,
        x[12] * x[36] +
        x[14] * x[39] +
        x[16] * x[42] +
        x[18] * x[45] +
        0.333333333333333 * x[20] == 30.0
    )
    @NLconstraint(
        m,
        e26,
        x[13] * x[34] +
        x[15] * x[37] +
        x[17] * x[40] +
        x[19] * x[43] +
        0.333333333333333 * x[21] == 70.0
    )
    @NLconstraint(
        m,
        e27,
        x[13] * x[35] +
        x[15] * x[38] +
        x[17] * x[41] +
        x[19] * x[44] +
        0.333333333333333 * x[21] == 50.0
    )
    @NLconstraint(
        m,
        e28,
        x[13] * x[36] +
        x[15] * x[39] +
        x[17] * x[42] +
        x[19] * x[45] +
        0.333333333333333 * x[21] == 70.0
    )
    @NLconstraint(m, e29, x[4] * x[28] - x[22] == 0.0)
    @NLconstraint(m, e30, x[4] * x[29] - x[23] == 0.0)
    @NLconstraint(m, e31, x[4] * x[30] - x[24] == 0.0)
    @NLconstraint(m, e32, x[5] * x[31] - x[25] == 0.0)
    @NLconstraint(m, e33, x[5] * x[32] - x[26] == 0.0)
    @NLconstraint(m, e34, x[5] * x[33] - x[27] == 0.0)
    @constraint(m, e35, x[34] + x[35] + x[36] == 1.0)
    @constraint(m, e36, x[37] + x[38] + x[39] == 1.0)
    @constraint(m, e37, x[40] + x[41] + x[42] == 1.0)
    @constraint(m, e38, x[43] + x[44] + x[45] == 1.0)
    @constraint(m, e39, x[28] + x[29] + x[30] == 1.0)
    @constraint(m, e40, x[31] + x[32] + x[33] == 1.0)
    @constraint(m, e41, x[36] == 0.0)
    @constraint(m, e42, x[43] == 0.0)
    @constraint(m, e43, x[4] - 300 * b[50] <= 0.0)
    @constraint(m, e44, x[5] - 300 * b[51] <= 0.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)
    return m
end

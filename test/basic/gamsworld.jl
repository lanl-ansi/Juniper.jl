function batch_problem()
    m = Model()
    @variable(m, 0 <= x1 <= 1.38629436111989)
    @variable(m, 0 <= x2 <= 1.38629436111989)
    @variable(m, 0 <= x3 <= 1.38629436111989)
    @variable(m, 0 <= x4 <= 1.38629436111989)
    @variable(m, 0 <= x5 <= 1.38629436111989)
    @variable(m, 0 <= x6 <= 1.38629436111989)
    @variable(m, 5.7037824746562 <= x7 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x8 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x9 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x10 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x11 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x12 <= 8.00636756765025)
    @variable(m, 4.45966 <= x13 <= 397.747)
    @variable(m, 3.7495 <= x14 <= 882.353)
    @variable(m, 4.49144 <= x15 <= 833.333)
    @variable(m, 3.14988 <= x16 <= 638.298)
    @variable(m, 3.04452 <= x17 <= 666.667)
    @variable(m, 0.729961 <= x18 <= 2.11626)
    @variable(m, 0.530628 <= x19 <= 1.91626)
    @variable(m, 1.09024 <= x20 <= 2.47654)
    @variable(m, -0.133531 <= x21 <= 1.25276)
    @variable(m, 0.0487901 <= x22 <= 1.43508)
    @variable(m, b23, Bin)
    @variable(m, b24, Bin)
    @variable(m, b25, Bin)
    @variable(m, b26, Bin)
    @variable(m, b27, Bin)
    @variable(m, b28, Bin)
    @variable(m, b29, Bin)
    @variable(m, b30, Bin)
    @variable(m, b31, Bin)
    @variable(m, b32, Bin)
    @variable(m, b33, Bin)
    @variable(m, b34, Bin)
    @variable(m, b35, Bin)
    @variable(m, b36, Bin)
    @variable(m, b37, Bin)
    @variable(m, b38, Bin)
    @variable(m, b39, Bin)
    @variable(m, b40, Bin)
    @variable(m, b41, Bin)
    @variable(m, b42, Bin)
    @variable(m, b43, Bin)
    @variable(m, b44, Bin)
    @variable(m, b45, Bin)
    @variable(m, b46, Bin)

    @NLobjective(
        m,
        Min,
        250 * exp(x1 + 0.6 * x7) +
        250 * exp(x2 + 0.6 * x8) +
        250 * exp(x3 + 0.6 * x9) +
        250 * exp(x4 + 0.6 * x10) +
        250 * exp(x5 + 0.6 * x11) +
        250 * exp(x6 + 0.6 * x12)
    )

    @NLconstraint(
        m,
        250 * exp(x1 + 0.6 * x7) +
        250 * exp(x2 + 0.6 * x8) +
        250 * exp(x3 + 0.6 * x9) +
        250 * exp(x4 + 0.6 * x10) +
        250 * exp(x5 + 0.6 * x11) +
        250 * exp(x6 + 0.6 * x12) <= 290000
    )

    @constraint(m, x7 - x13 >= 2.06686275947298)
    @constraint(m, x8 - x13 >= 0.693147180559945)
    @constraint(m, x9 - x13 >= 1.64865862558738)
    @constraint(m, x10 - x13 >= 1.58923520511658)
    @constraint(m, x11 - x13 >= 1.80828877117927)
    @constraint(m, x12 - x13 >= 1.43508452528932)
    @constraint(m, x7 - x14 >= -0.356674943938732)
    @constraint(m, x8 - x14 >= -0.22314355131421)
    @constraint(m, x9 - x14 >= -0.105360515657826)
    @constraint(m, x10 - x14 >= 1.22377543162212)
    @constraint(m, x11 - x14 >= 0.741937344729377)
    @constraint(m, x12 - x14 >= 0.916290731874155)
    @constraint(m, x7 - x15 >= -0.356674943938732)
    @constraint(m, x8 - x15 >= 0.955511445027436)
    @constraint(m, x9 - x15 >= 0.470003629245736)
    @constraint(m, x10 - x15 >= 1.28093384546206)
    @constraint(m, x11 - x15 >= 1.16315080980568)
    @constraint(m, x12 - x15 >= 1.06471073699243)
    @constraint(m, x7 - x16 >= 1.54756250871601)
    @constraint(m, x8 - x16 >= 0.832909122935104)
    @constraint(m, x9 - x16 >= 0.470003629245736)
    @constraint(m, x10 - x16 >= 0.993251773010283)
    @constraint(m, x11 - x16 >= 0.182321556793955)
    @constraint(m, x12 - x16 >= 0.916290731874155)
    @constraint(m, x7 - x17 >= 0.182321556793955)
    @constraint(m, x8 - x17 >= 1.28093384546206)
    @constraint(m, x9 - x17 >= 0.8754687373539)
    @constraint(m, x10 - x17 >= 1.50407739677627)
    @constraint(m, x11 - x17 >= 0.470003629245736)
    @constraint(m, x12 - x17 >= 0.741937344729377)
    @constraint(m, x1 + x18 >= 1.85629799036563)
    @constraint(m, x2 + x18 >= 1.54756250871601)
    @constraint(m, x3 + x18 >= 2.11625551480255)
    @constraint(m, x4 + x18 >= 1.3609765531356)
    @constraint(m, x5 + x18 >= 0.741937344729377)
    @constraint(m, x6 + x18 >= 0.182321556793955)
    @constraint(m, x1 + x19 >= 1.91692261218206)
    @constraint(m, x2 + x19 >= 1.85629799036563)
    @constraint(m, x3 + x19 >= 1.87180217690159)
    @constraint(m, x4 + x19 >= 1.48160454092422)
    @constraint(m, x5 + x19 >= 0.832909122935104)
    @constraint(m, x6 + x19 >= 1.16315080980568)
    @constraint(m, x1 + x20 >= 0)
    @constraint(m, x2 + x20 >= 1.84054963339749)
    @constraint(m, x3 + x20 >= 1.68639895357023)
    @constraint(m, x4 + x20 >= 2.47653840011748)
    @constraint(m, x5 + x20 >= 1.7404661748405)
    @constraint(m, x6 + x20 >= 1.82454929205105)
    @constraint(m, x1 + x21 >= 1.16315080980568)
    @constraint(m, x2 + x21 >= 1.09861228866811)
    @constraint(m, x3 + x21 >= 1.25276296849537)
    @constraint(m, x4 + x21 >= 1.19392246847243)
    @constraint(m, x5 + x21 >= 1.02961941718116)
    @constraint(m, x6 + x21 >= 1.22377543162212)
    @constraint(m, x1 + x22 >= 0.741937344729377)
    @constraint(m, x2 + x22 >= 0.916290731874155)
    @constraint(m, x3 + x22 >= 1.43508452528932)
    @constraint(m, x4 + x22 >= 1.28093384546206)
    @constraint(m, x5 + x22 >= 1.30833281965018)
    @constraint(m, x6 + x22 >= 0.78845736036427)

    @NLconstraint(
        m,
        250000 * exp(x18 - x13) +
        150000 * exp(x19 - x14) +
        180000 * exp(x20 - x15) +
        160000 * exp(x21 - x16) +
        120000 * exp(x22 - x17) <= 6000
    )

    @NLconstraint(
        m,
        x1 - 0.693147180559945 * b29 - 1.09861228866811 * b35 -
        1.38629436111989 * b41 == 0
    )

    @NLconstraint(
        m,
        x2 - 0.693147180559945 * b30 - 1.09861228866811 * b36 -
        1.38629436111989 * b42 == 0
    )

    @NLconstraint(
        m,
        x3 - 0.693147180559945 * b31 - 1.09861228866811 * b37 -
        1.38629436111989 * b43 == 0
    )

    @NLconstraint(
        m,
        x4 - 0.693147180559945 * b32 - 1.09861228866811 * b38 -
        1.38629436111989 * b44 == 0
    )

    @NLconstraint(
        m,
        x5 - 0.693147180559945 * b33 - 1.09861228866811 * b39 -
        1.38629436111989 * b45 == 0
    )

    @NLconstraint(
        m,
        x6 - 0.693147180559945 * b34 - 1.09861228866811 * b40 -
        1.38629436111989 * b46 == 0
    )

    @constraint(m, b23 + b29 + b35 + b41 == 1)
    @constraint(m, b24 + b30 + b36 + b42 == 1)
    @constraint(m, b25 + b31 + b37 + b43 == 1)
    @constraint(m, b26 + b32 + b38 + b44 == 1)
    @constraint(m, b27 + b33 + b39 + b45 == 1)
    @constraint(m, b28 + b34 + b40 + b46 == 1)
    return m
end

function cvxnonsep_nsig20r_problem()
    m = Model()
    @variable(m, 1 <= i1 <= 10, Int)
    @variable(m, 1 <= i2 <= 10, Int)
    @variable(m, 1 <= i3 <= 10, Int)
    @variable(m, 1 <= i4 <= 10, Int)
    @variable(m, 1 <= i5 <= 10, Int)
    @variable(m, 1 <= i6 <= 10, Int)
    @variable(m, 1 <= i7 <= 10, Int)
    @variable(m, 1 <= i8 <= 10, Int)
    @variable(m, 1 <= i9 <= 10, Int)
    @variable(m, 1 <= i10 <= 10, Int)
    @variable(m, 1E-5 <= x11 <= 10)
    @variable(m, 1E-5 <= x12 <= 10)
    @variable(m, 1E-5 <= x13 <= 10)
    @variable(m, 1E-5 <= x14 <= 10)
    @variable(m, 1E-5 <= x15 <= 10)
    @variable(m, 1E-5 <= x16 <= 10)
    @variable(m, 1E-5 <= x17 <= 10)
    @variable(m, 1E-5 <= x18 <= 10)
    @variable(m, 1E-5 <= x19 <= 10)
    @variable(m, 1E-5 <= x20 <= 10)
    @variable(m, x21)
    @variable(m, x22)
    @variable(m, x23)
    @variable(m, x24)
    @variable(m, x25)
    @variable(m, x26)
    @variable(m, x27)
    @variable(m, x28)
    @variable(m, x29)
    @variable(m, x30)
    @variable(m, x31)
    @variable(m, x32)
    @variable(m, x33)
    @variable(m, x34)
    @variable(m, x35)
    @variable(m, x36)
    @variable(m, x37)
    @variable(m, x38)
    @variable(m, x39)
    @variable(m, x40)
    @constraint(
        m,
        x21 +
        x22 +
        x23 +
        x24 +
        x25 +
        x26 +
        x27 +
        x28 +
        x29 +
        x30 +
        x31 +
        x32 +
        x33 +
        x34 +
        x35 +
        x36 +
        x37 +
        x38 +
        x39 +
        x40 <= -1.6094379124341
    )
    @NLconstraint(m, -0.065 * log(i1) - x21 <= 0)
    @NLconstraint(m, -0.004 * log(i2) - x22 <= 0)
    @NLconstraint(m, -0.084 * log(i3) - x23 <= 0)
    @NLconstraint(m, -0.093 * log(i4) - x24 <= 0)
    @NLconstraint(m, -0.06 * log(i5) - x25 <= 0)
    @NLconstraint(m, -0.075 * log(i6) - x26 <= 0)
    @NLconstraint(m, -0.074 * log(i7) - x27 <= 0)
    @NLconstraint(m, -0.039 * log(i8) - x28 <= 0)
    @NLconstraint(m, -0.065 * log(i9) - x29 <= 0)
    @NLconstraint(m, -0.017 * log(i10) - x30 <= 0)
    @NLconstraint(m, -0.07 * log(x11) - x31 <= 0)
    @NLconstraint(m, -0.03 * log(x12) - x32 <= 0)
    @NLconstraint(m, -0.028 * log(x13) - x33 <= 0)
    @NLconstraint(m, -0.005 * log(x14) - x34 <= 0)
    @NLconstraint(m, -0.01 * log(x15) - x35 <= 0)
    @NLconstraint(m, -0.082 * log(x16) - x36 <= 0)
    @NLconstraint(m, -0.069 * log(x17) - x37 <= 0)
    @NLconstraint(m, -0.032 * log(x18) - x38 <= 0)
    @NLconstraint(m, -0.095 * log(x19) - x39 <= 0)
    @NLconstraint(m, -0.003 * log(x20) - x40 <= 0)
    @NLobjective(
        m,
        Min,
        1.5 * i1 +
        0.51 * i2 +
        1.01 * i3 +
        1.4 * i4 +
        1.78 * i5 +
        1.92 * i6 +
        1.09 * i7 +
        0.48 * i8 +
        0.67 * i9 +
        0.52 * i10 +
        1.68 * x11 +
        0.51 * x12 +
        1.63 * x13 +
        0.49 * x14 +
        1.86 * x15 +
        0.7 * x16 +
        0.39 * x17 +
        0.5 * x18 +
        1.23 * x19 +
        0.95 * x20
    )

    return m
end

using Distances

function load_fac(file_name)
    f = open("./test/data/" * file_name)
    lines = readlines(f)
    n_str,m_str = split(lines[1])
    N, M = parse(Int, n_str), parse(Int, m_str)

    f_s = Vector{Float64}(N)
    f_c = Vector{Float64}(N)
    f_pos = [Vector{Float64}(2) for _ in 1:N]

    for i = 2:N+1
        s_str, c_str, x_str, y_str = split(lines[i])
        f_s[i-1], f_c[i-1] = parse(Float64, s_str), parse(Float64, c_str)
        f_pos[i-1] = [parse(Float64, x_str), parse(Float64, y_str)]
    end

    f_c_max,_ = findmax(f_c)
    f_c /= f_c_max

    f_s_max,_ = findmax(f_s)
    f_s /= f_s_max

    c_d = Vector{Float64}(M)
    c_pos = [Vector{Float64}(2) for _ in 1:M]

    for i = N+2:length(lines)
        d_str, x_str, y_str = split(lines[i])
        c_d[i-N-1] = parse(Float64, d_str)
        c_pos[i-N-1] = [parse(Float64, x_str), parse(Float64, y_str)]
    end

    c_d /= f_c_max

    dist_mat = zeros(Float64, M, N)
    for i=1:N, j=1:M
        dist_mat[j, i] = euclidean(c_pos[j],f_pos[i])
    end

    dist_mat_max,_ = findmax(dist_mat)
    dist_mat /= dist_mat_max

    return N,M,f_s,f_c,f_pos,c_d,c_pos,dist_mat
end
using Oscar
using TensorCast

lyndonlist23 = [0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 7, 8, 11, 13, 17,  21, 28, 34, 45,
56, 73, 92, 120, 151, 197, 250, 324, 414, 537, 687, 892, 1145, 1484, 1911, 2479,
3196, 4148, 5359, 6954, 9000, 11687, 15140, 19672, 25516]

lyndonlist2301 = [0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

lyndonlistNN = [1,1,2,3,6,9,18, 30, 56, 99, 186,335]
lyndonlistNN0 = [0,1,2,3,6,9,18, 30, 56, 99, 186,335]
lyndonlistNN001 = [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
#https://oeis.org/A059966
#1,2,12,3,112,13,4,1112,122,113,14,23,5,11112,1122,1113,114,123,132,15,24,6

#1,2,4,8
#1,11,2,111,12,21,3,1111,112,121,211,22,13,31,4,1112,1121,1211,2111,122,212,221,
function combinat(li,k)
    if length(li)==1
        in = k/li[1]
        if isinteger(in) == true
            return [Int(in)]
        else
            return []
        end
    end
    combs=[]
    l = length(li)
    for red in 0:Int(floor(k/li[l]))
        for c in combinat(li[1:end-1],k-red*li[l])
            combs = [combs;[[c;[red]]]]
        end
    end
    return combs
end

function combinat_mzv_var23(li,k)
    lis=li[1:k]
    #print(lis)
    newli =[]
    j=1
    for i in lis
        newli=[newli;[j for l in 1:i]]
        j=j+1
    end
    #print(newli)
    return convert(Vector{Vector{Int64}},combinat(newli,k))
end

function red_combinat_mzv_var23_old(li,k)
    C = combinat_mzv_var23(li,k)
    C = C[1:end-1]
    red = Int64[]
    for i in 1:length(C[1])
        for j in 1:length(C)
            if C[j][i] != 0
                red = push!(red,i)
                break
            end
        end
    end
    return [c[red] for c in C]
end

function red_combinat_mzv_var23(li,k)
    C = combinat_mzv_var23(li,k)
    red = Int64[]
    for i in 1:length(C[1])
        for j in 1:length(C)
            if C[j][i] != 0
                red = push!(red,i)
                break
            end
        end
    end
    return [c[red] for c in C]
end

function vred_combinat_mzv_var23(li,k)
    C = red_combinat_mzv_var23(li,k)
    del1 = Int64[]
    del2 = Int64[]
    k = Int64(0)
    n = Int64(0)
    for i in 1:length(C[1])
        for j in 1:length(C)
            if C[j][i] != 0
                n = n+1
                k = j
            end
        end
        if n == 1
            push!(del2,i)
            push!(del1,k)
        end
        n = 0
    end
    red1 = Int64[]
    red2 = Int64[]
    for i in 1:length(C[1])
        if !(i in del2)
            push!(red2,i)
        end
    end
    for j in 1:length(C)
        if !(j in del1)
            push!(red1,j)
        end
    end
    return [c[red2] for c in C[red1]]
end

function kernelcheck_vred_combinat_mzv_var23(li,k)
    K = kernel(matrix(QQ,vred_combinat_mzv_var23(li,k)))
    for j in 1:size(K,2)
        if iszero(K[:,j]) == true
            return false
        end
    end
    return true
end


function lattice_mzv_var23(li,k)
    C = red_combinat_mzv_var2old(li,k)
    L = integer_lattice(matrix(ZZ,[C[1]]))
    for c in C[2:end]
        L = L + integer_lattice(matrix(ZZ,[c]))
    end
    return L
end

function index_lattice_mzv_var23(li,k)
    L = lattice_mzv_var23(li,k)
    K = integer_lattice(identity_matrix(ZZ,rank(L)))
    return index(K,L)
end

function lattice_vred_mzv_var23(li,k)
    C = vred_combinat_mzv_var23(li,k)
    Cdehom = map(x->x-C[end],C[1:end-1])
    Cproj = map(x->x[1:end-1],Cdehom)
    L = integer_lattice(matrix(ZZ,[Cproj[1]]))
    for c in Cproj[2:end]
        L = L + integer_lattice(matrix(ZZ,[c]))
    end
    return L
end

function index_vred_lattice_mzv_var23(li,k)
    L = lattice_vred_mzv_var23(li,k)
    K = integer_lattice(identity_matrix(ZZ,rank(L)))
    return index(K,L)
end

function red_polytope_mzv_var23(li,k)
    C = red_combinat_mzv_var23_old(li,k)
    C = push!(C,[Int64(0) for c in C[1]])
    return convex_hull(C)
end

function vred_polytope_mzv_var23(li,k)
    C = vred_combinat_mzv_var23(li,k)
    return convex_hull(C)
end

function vred_polytope_mzv_var23_in_lattice(li,k)
    C = vred_combinat_mzv_var23(li,k)
    Cdehom = map(x->x-C[end],C)
    Cproj = map(x->x[1:end-1],Cdehom)
    L = lattice_vred_mzv_var23(li,k)
    M = transpose(inv(basis_matrix(L)))
    Cp = map(x->M*x,Cproj)
    return convex_hull(Cp)
end

function matrix_mzv_var23(li,k)
    lis=li[1:k]
    #print(lis)
    newli =[]
    j=1
    for i in lis
        newli=[newli;[j for l in 1:i]]
        j=j+1
    end
    #print(newli)
    return matrix(ZZ,convert(Vector{Vector{Int64}},combinat(newli,k)))
end

function matrix_vred_mzv_var23(li,k)
    return matrix(ZZ,convert(Vector{Vector{Int64}},vred_combinat_mzv_var23(li,k)))
end

function polytope_mzv_var23(li,k)
    lis=li[1:k]
    #print(lis)
    newli =[]
    j=1
    for i in lis
        newli=[newli;[j for l in 1:i]]
        j=j+1
    end
    #print(newli)
    return convex_hull(convert(Vector{Vector{Int64}},combinat(newli,k)))
end

function deg_mzv_var23(li,k)
    lis=li[1:k]
    #print(lis)
    newli =[]
    j=1
    for i in lis
        newli=[newli;[j for l in 1:i]]
        j=j+1
    end
    #print(newli)
    print(convert(Vector{Vector{Int64}},combinat(newli,k)))
    return lattice_volume(convex_hull(convert(Vector{Vector{Int64}},combinat(newli,k))))
end

function vol_simplex(li,k)
    lis=li[1:k]
    #print(lis)
    newli =[]
    j=1
    for i in lis
        newli=[newli;[j for l in 1:i]]
        j=j+1
    end
    vertices = []
    j=1
    N = length(newli)
    print(newli)
    for i in newli
        newv = [QQ(0) for l in 1:N]
        newv[j]= QQ(k//i)
        vertices=[vertices;[newv]]
        j=j+1
    end
    vertices=[vertices;[[QQ(0) for i in 1:N]]]
    return convex_hull(convert(Vector{Vector{QQFieldElem}},vertices))
end

function check_smoothness(li,k)
    v=[collect(i) for i in vertices(convex_hull(vred_combinat_mzv_var23(li,k)))]
    for i in 1:length(v)
        C=cone(map(x->x-v[i],v[[1:i-1;i+1:length(v)]]))
        print(length(rays(C))," ")
    end
end

function degree_formula(li,k)
    return volume(red_polytope_mzv_var23(li,k))*factorial(big(dim(polytope_mzv_var23(li,k))))/index_lattice_mzv_var23(li,k)
end

function toric_ideal_mzv_var23(li,k)
    M = matrix_mzv_var23(li,k)
    n = number_of_rows(M)
    vars = ["x"*string(i) for i in 1:n]
    R,_ = graded_polynomial_ring(QQ,vars,[1 for i in 1:n])
    return toric_ideal(R,M)
end

function toric_ideal_vred_mzv_var23(li,k)
    M = matrix_vred_mzv_var23(li,k)
    n = number_of_rows(M)
    vars = ["x"*string(i) for i in 1:n]
    R,_ = graded_polynomial_ring(QQ,vars,[1 for i in 1:n])
    return toric_ideal(R,M)
end

function int_deg(v)
    return map(x->parse(Int64,repr(degree(x))[begin+1:end-1]),v)
end

#vertices(convex_hull(lattice_points(vol_simplex(lyndonlist23,18))))

# M1= transpose([6 3 2 2 1 1 0 0; 0 2 1 0 1 0 1 0; 0 0 1 0 0 0 0 0; 0 0 0 0 1 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 0 0 1 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1]
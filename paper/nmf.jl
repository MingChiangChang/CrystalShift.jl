# Non-negative Matrix Factorization via Fast Conical Hull Algorithm
# Adapted from https://github.com/SebastianAment/PhaseMapping.jl/blob/master/src/nmf.jl
using LinearAlgebra

############################# helper functions #################################
colnorms(X::AbstractMatrix) = [norm(x) for x in eachcol(X)]
# matrix-matrix multiply and projection onto positive values
function posmul!(AB, A, B)
    mul!(AB, A, B)
    return @. AB = max(AB, 0)
end

################################################################################
# Structure holding temporaries for XRAY algorithm
struct XRAY{T, XT<:AbstractMatrix{T}, KT<:AbstractVector{Int}, ST}
    X::XT # data matrix
    R::XT # residual matrix
    H::XT # activation matrix
    K::KT
    selection::ST # ray selection algorithm
end

# 1: dist ray, 2: greedy, 3: max, 4: rand
function XRAY(X::AbstractMatrix, r::Int, selection_id::Int = 3)
    begin # ray selection strategies
        colnormsX = colnorms(X)
        m = size(X, 2)
        RX = similar(X, (m, m)) # temporary for R'X
        dist(R, X) = argmax(colnorms(posmul!(RX, R', X)))
        greedy(R, X) = argmax(colnorms(posmul!(RX, R', X)) ./ colnormsX)
        maxres(R, X) = argmax(colnorms(R)) # based on maximum residual
        randres(R, X) = rand(findall(>(0), colnorms(R)))
        selection_algorithms = (dist, greedy, maxres, randres)
    end
    selection = selection_algorithms[selection_id]
    XRAY(X, r, selection)
end
function XRAY(X::AbstractMatrix, r::Int, selection::Function)
    R = copy(X)
    H = zeros(eltype(X), (r, size(X, 2)))
    K = zeros(Int, 0) # anchor index set
    XRAY(X, R, H, K, selection)
end

################################ xray algorithm ################################
function xray(X::AbstractMatrix, r::Int, tol::Real = 1e-6; selection_id::Int = 3)
    U = XRAY(X, r, selection_id)
    return xray!(U, r, tol)
end

# identify anchor columns of X
# r is maximum rank of NMF
function xray!(U::XRAY, r::Int, tol::Real = 1e-6)
    X, R, H, K = U.X, U.R, U.H, U.K
    n, m = size(X)
    i = 0
    while norm(R) > tol && i < r
        i += 1
        k = U.selection(R, X)

        # extend index set
        push!(K, k)
        Xk, Hi = @views X[:, K], H[1:i, :]

        # projection, solves min | X - X[:,I]*H |_2 s.t. H ≥ 0
        snnls!(Hi, Xk, X) # nnls!(Hi, Xk, X)

        # update residual
        @. R = X
        @views mul!(R, Xk, Hi, -1, 1)
    end
    W, H = @views X[:, K], H[1:i, :]
    return W, H, K
end

# using LinearAlgebraExtensions: LowRank
# # alternating non-negative least squares NMF
# could be used after XRAY for refinement
# function annls!(L::LowRank, X::AbstractMatrix; tol::Real = 1e-6, maxiter::Int = 128)
#     W, H = L.U, L.V
#     R = copy(X)
#     mul!(R, W, H, -1, 1)
#     obj = norm(R)
#     for i in 1:maxiter
#         snnls!(H, W, X) # solve min_H |X-WH| s.t. H ≥ 0
#         snnls!(W', H', X') # solve min_W' |X'-H'W'| s.t. W' ≥ 0
#         mul!(R, W, H, -1, 1)
#         newobj = norm(R)
#         obj-newobj > tol || break
#         obj = newobj
#     end
#     return L
# end

############################### NNLS solvers ###################################
# simultaneous nnls: min |AX - B|
# TODO: pre-allocate all temporaries in XRAY
function snnls!(X::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix;
            maxiter::Int = 128, tol::Real = 1e-6)
    m, n = size(B)
    r = size(A, 2)
    BA = B'A
    AA = A'A
    XAA = X'*AA
    R = similar(B)
    obj = Inf
    xi = zeros(n)
    for j in 1:maxiter
        for i in 1:size(X, 1)
            Xi, XAAi, BAi, AAi = @views X[i,:], XAA[:,i], BA[:,i], AA[:,i]
            copy!(xi, Xi) # store old value for update
            @. Xi = max(Xi - (XAAi - BAi) / AA[i,i], 0)
            @. XAA += (Xi-xi) * AAi' # outer product
        end
        @. R = B
        mul!(R, A, X, -1, 1) # update residual
        newobj = norm(R)
        obj-newobj > tol || break
        obj = newobj
    end
    return X
end

# slower and adds dependencies, though potentially better for high accuracy
# using JuMP
# using JuMP: SecondOrderCone
# using ECOS # ECOS is interior point method (more accurate but slower on large problems)
# # using SCS # operator splitting method, fast but not high accuracy
# function nnls(A::AbstractMatrix, b::AbstractVector)
#     model = JuMP.Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => false))
#     n, m = size(A)
#     @variable(model, x[1:m] ≥ 0)
#     @variable(model, r[1:n])
#     @constraint(model, ctr[i in 1:n], r[i] == b[i]-sum(A[i,j]*x[j] for j in 1:m)) # r = y-A*x
#     @objective(model, Min, sum(r[i]*r[i] for i in 1:n))
#     # return model
#     JuMP.optimize!(model)
#     return JuMP.value.(x)
# end
#
# function nnls(A::AbstractMatrix, B::AbstractMatrix)
#     X = zeros(size(A, 2), size(B, 2))
#     return nnls!(X, A, B)
# end
# function nnls!(X::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
#     for i in 1:size(B, 2)
#         X[:, i] = nnls(A, @view B[:, i])
#     end
#     return X
# end
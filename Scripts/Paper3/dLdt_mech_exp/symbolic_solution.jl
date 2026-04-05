using Symbolics

M = 2

function tri_matrix_dense_free(N::Int)
    A = zeros(Float64, N, N)
    for i in 1:N
        A[i,i] = -2
        if i < N
            A[i, i+1] = 1
            A[i+1, i] = 1
        end
    end
    A[N, N] = -1 # free boundary condition modification
    return A
end


@variables x1(t) x2(t) x1_0 x2_0 t k η a
#k_c = 1.0; η_c = 1.0; a_c = 2.0;
A = tri_matrix_dense_free(M)
DiagA = eigen(A)
λ = DiagA.values
S = DiagA.vectors
#S = [1.0 -2.0 - λ[1]; 1.0 -2.0 - λ[2]]
S_inv = inv(S)
E = Matrix{Num}(undef, M, M)
for i in 1:M
    for j in 1:M
        E[i,j] = 0.0
    end
    E[i,i] = exp((k/η * λ[i] * t))
end
y_vec = zeros(Num,M); y_vec[end] = a;
z0 =  x0 + inv(A) * y_vec
x = S * E * S_inv * z0
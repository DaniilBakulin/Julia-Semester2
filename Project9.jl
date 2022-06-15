function upper_triangular(N::Integer)
   A = zeros(N, N)
   for i in 1:N
        A[i, i:N] = N * randn(1, N-i+1)
        if A[i, i]==0
            A[i, i] = 1.0
        end
    end
    return A
end

function reverse_gauss_1(A::AbstractMatrix{Float64}, b::AbstractVector{Float64})
    x = similar(b)
    N = size(A, 1)
    for k in 0:N-1
        sm = 0
        x[N-k]=(b[N-k] - dot(@view(A[N-k, N-k+1:N]), @view(x[N-k+1:N])))/A[N-k, N-k]

    end
    return x
end



function transform_to_steps!(A::AbstractMatrix; epsilon = 1e-7, degenerate_exception = true)
    for k ∈ 1:size(A, 1)
         absval, Δk = findmax(abs, @view(A[k:end,k]))
         (degenerate_exception && absval <= epsilon) && throw("Вырожденая матрица")
         if Δk > 1
             A[k,k:end], A[k+Δk-1,k:end] = A[k+Δk-1,k:end], A[k,k:end]
         end 
         for i ∈ k+1:size(A,1)
             t = A[i,k]/A[k,k]          
             @. @views A[i,k:end] = A[i,k:end] - t * A[k,k:end] 
         end
     end
 end
using LinearAlgebra

#Task_1
function transform_to_steps!(A::AbstractMatrix{Float64}, epsilon = 1e-7, degenerate_exception = true)
    for k ∈ 1:size(A, 1)
         absval, Δk = findmax(abs, @view(A[k:end,k]))
         (degenerate_exception && absval <= epsilon) && throw("Вырожденая матрица")
         Δk > 1 && swap!(@view(A[k,k:end]), @view(A[k+Δk-1,k:end])) 
         for i ∈ k+1:size(A,1)
             t = A[i,k]/A[k,k]          
             @. @views A[i,k:end] = A[i,k:end] - t * A[k,k:end] 
         end
     end
     return A
 end
 
 function swap!(A,B)
    for i in eachindex(A)
        A[i], B[i] = B[i], A[i]
    end
end

#Task_2
function transform_to_steps2!(A::AbstractMatrix{Float64}, epsilon = 1e-7, degenerate_exception = true)
    for k ∈ size(A, 1):-1:1
         absval, Δk = findmax(abs, @view(A[k,1:k]))
         (degenerate_exception && absval <= epsilon) && throw("Вырожденая матрица")
         Δk > 1 && swap!(@view(A[k:-1:1,k]), @view(A[k:-1:1,k-Δk+1])) 
         for i ∈ k-1:-1:1
             t = A[k,i]/A[k,k]          
             @. @views A[k:-1:1,i] = A[k:-1:1,i] - t * A[k:-1:1,k] 
         end
     end
     return A
end 
#Task_3

function diag_zero_element(A::AbstractMatrix{Float64}, epsilon = 1e-7, degenerate_exception = false)
    C = []
    AA =deepcopy(A)
    transform_to_steps!(AA)
    for k in 1::size(A, 1)
        if (abs(AA[k, k])<=epsilon  )
            push!(C,k)
        end
    end
    return C
end
#Task_4

function rang_(A::AbstractMatrix)
    B = diag_zero_element(A, 1e-7,false)
    return size(A, 1) - size(B, 1)
end

function det_(A::AbstractMatrix)#fixed
    AA=deepcopy(A)
    k = 1
    transform_to_steps!(AA,1e-7)
    for i in 1:size(A, 1)
        k=k*AA[i, i]
    end
    return k
end

function detp(A::AbstractMatrix)
    sm = 0
    x = 0
    if size(A, 1)==1
        return A[1, 1]
    else
        for x in 1:size(A, 1)
            x = A[1, x] * ((-1)^(1+x)) * detp(slice(A, 1, x))
            if x != nothing
                sm+=x
            end
        end
        return sm
    end 
end
#Task_5
function reverse_gauss_1(A::AbstractMatrix{Float64}, b::AbstractVector{Float64})
    x = similar(b)
    N = size(A, 1)
    for k in 0:N-1
        sm = 0
        x[N-k]=(b[N-k] - dot(@view(A[N-k, N-k+1:N]), @view(x[N-k+1:N])))/A[N-k, N-k]

    end
    return x
end

solve_Axb(A, b) = reverse_gauss_1(transform_to_steps!(deepcopy(A)), deepcopy(b))
#Task_6
#Task_7

function slice(A::AbstractMatrix{Float64}, i, j)
    n = size(A, 1)
    B = zeros(n-1, n-1)
    i_ = 0
    j_ = 0
    for x in 1:n-1
        for y in 1:n-1
            i_ = x
            j_ = y
            if x >=i
                i_+=1
            end
            if y>=j
                j_+=1
            end
            B[x,y] = A[i_, j_]
        end
    end
    return B
end
function inv_(A::AbstractMatrix{Float64})
    n=size(A,1)
    AA = deepcopy(A)
    dt = detp(AA)
    for i in 1:n
        for j in 1:n
            AA[j, i] = ((-1)^(i+j))*detp(slice(A, i, j))
        end
    end
    return 1/dt*AA
end


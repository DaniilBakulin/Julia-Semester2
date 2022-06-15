function findall_(f::Function, A)
    s = []
    for i in 1:length(A)
        if f(A[i])
            push!(s, i)
        end
    end
    return s
end

function findfirst_(A, a)
    for i in 1:length(A)
        if A[i] == a
            return i 
        end
    end
    return nothing
end

function findlast_(A, a)
    ix = nothing
    for i in 1:length(A)
        if A[i] == a
            ix=i 
        end
    end
    return ix
end

function filter_(f, A)
    s = []
    for i in 1:length(A)
        if f(A[i]) != false
            push!(s, A[i])
        end
    end
    return s
end


function issorted(A)  
    for i in 1:length(A)-1
        if A[i+1]<A[i]
            return false
        end
    end
    return true
end

function f(x)
    return x
end

function bubblesort!(a, by::Function = f(x))
    n = length(a)
    for k in 1:n-1
        is_sorted = true
        for i in 1:n-k
            if by(a[i])>by(a[i+1])
                a[i], a[i+1] = a[i+1], a[i]
                is_sorted = false
            end
        end
        if is_sorted
            break
        end
    end
    return a
end
bubblesort(a) = bubblesort!(deepcopy(a))

function bubblesortperm!(a, by::Function = f(x))
    n = length(a)
    indexes = collect(1:n)
    for k in 1:n-1
        is_sorted = true
        for i in 1:n-k
            if by(a[i]) > by(a[i+1])
                a[i], a[i+1] = a[i+1], a[i]
                indexes[i], indexes[i+1] = indexes[i+1], indexes[i]
                is_sorted = false
            end
        end
        if is_sorted
            break
        end
    end
    return indexes
end
bubblesortperm(a) = bubblesortperm!(deepcopy(a))


function slice(A::Matrix,I::Vector{Int},J::Vector{Int})
    B=Matrix{eltype(A)}(undef,length(I),length(J))
    for i in I
        for j in J
            B[i,j]=A[I[i],J[j]]
        end
    end
    return B
end

function col_bubblesort!(A::AbstractMatrix)
    for j in size(A,2)
        bubblesort!(@view A[:,j])
    end
    return A
end

col_sortkey!(A::AbstractMatrix, key_values) = A[:, sortperm!(key_values)]

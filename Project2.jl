#Task_1
function issorted_(a)
    prev = a[1]
    for i in a[2:end]
        if i < prev
            return false
        end
        prev = i
    end
    return true
end

function binary_search(a, A, b, e)
    if b > e
        return false
    end
    i = (b+e+1) ÷ 2 

    if A[i] == a 
        return true
    end
    if a>A[i]
        return binary_search(a, A, i+1, e)
    else
        return binary_search(a, A, b, i-1)
    end
end

function length_(A)
    len = 0
    for _ in A
        len += 1
    end
    return len
end

function sum_(A::AbstractVector{T}) where T
    s = T(0)
    for a in A
        s += a
    end
    return s
end


function prod_(A)
    p = eltype(A)(1)
    for a in A
        p *= a
    end
    return p
end

function maximum_(A)
    M = typemin(eltype(A)) # m = -Inf
    for a in A
        M = max(M,a)
    end
    return M
end

function minimum_(A)
    m = typemin(eltype(A)) # m = -Inf
    for a in A
        m = min(m,a)
    end
    return m
end

function imax_(A)
    @assert !isempty(A)
    imax = firstindex(A)
    for k in eachindex(A)
        if A[k] > A[imax] 
            imax = k
        end
    end
    return imax
end

function insertsort!(A)
    n=length(A)
    for k in eachindex(A) #2:n
        # часть массива A[1:k-1] уже отсортирована
        op_insert!(A,k)
    end
    return A
end

op_insert!(A,k) =
    while k>1 && A[k-1] > A[k]
        A[k-1], A[k] = A[k], A[k-1]
        k -= 1
    end

function evalpoly_(x,A)
    Q = first(A) # коэффициенты записаны по убыванию степени
    for a in @view A[2:end]
        Q=Q*x+a
    end
    return Q
end

function mean_(A::AbstractVector{T}) where T
    s = T(0)
    i = 0
    for a in A
        i += 1
        s += a
    end
    return s / i
end

function count_maximum_(A)
    M = typemin(eltype(A)) # m = -Inf
    c = 0
    for a in A
        if a == M
            c+=1
        end
        if a > M
            M = a
            c = 1
        end
    end
    return (c, M)
end

function _isperm(p)
    n = length(p)
    used = falses(n) # возвращает нулевой BitVector длины n
    for i in p
        if i<0 || i>length(p)
            return false
        elseif used[i]
            return false
        end
        used[i] = true 
    end
    return true
end

function _invpermute!(A, p)
    for i in p
        if i > 0
            A[i], A[p[i]] = A[p[i]], A[i]
            p[i] = -p[i]
        end
    end
    for i in eachindex(p)
        p[i] = -p[i]
    end
    return A   
end

function _permute!(A, p) 
    for i in eachindex(p)
        if p[i] < 0
            continue
        end 
        buff = A[i]
        j_prev, j = i, p[i]           
        p[i] = -p[i]
        while j != i                
            A[j_prev] = A[j]
            j_prev, j = j, p[j]            
            p[j_prev] = -p[j_prev]
        end        
        A[j_prev] = buff 
        # перемещения элементов массива A по очередному циклу (по очередной циклической перестановке индексов) полностью завершены
    end
    for i in eachindex(p)
        p[i] = -p[i]
    end        
    return A
end

#Task_2
function _reverse!(a)
    for i in 1:length(a) ÷ 2
        a[i], a[end-i+1] = a[end-i+1], a[i] 
    end
end
#Task_3
function shift_right!(a)
    resize!(a, length(a)+1)
    a[begin+1:end] = @view(a[begin:end-1])
    a[begin] = 0
end

#Task_4
function circshift_right!(a)
    buf = a[end]
    for i in lastindex(a):-1:firstindex(a)-1
        a[i] = a[i-1]
    end
    a[begin] = buf
end

function circshift_left!(a)
    buf = a[begin]
    for i in firstindex(a):lastindex(a)-1
        a[i] = a[i+1]
    end
    a[begin] = buf
end

#Task_5
function _circshift!(a, k)
    if k > 0
        reverse!(@view(a[1:k]))
        reverse!(@view(a[k+1:end]))
        reverse!(a) 
    elseif k < 0
        reverse!(a)
        reverse!(@view(a[1:-k]))
        reverse!(@view(a[-k+1:end]))
    end
    return a
end

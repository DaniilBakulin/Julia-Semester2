#Task_1
function quick_sort!(A)
    if isempty(A)
        return A
    end
    N = length(A)
    K, M = part_sort!(A, A[rand(1:N)]) # - "базовый" элемент массива выбирается случайным образом
    quick_sort!(@view A[1:K])
    quick_sort!(@view A[M:N])
    return A
end

function part_sort!(A, b)
    N = length(A)
    K, L, M = 0, 0, N
    @inbounds while L < M
        if A[L+1] == b
            L += 1
        elseif A[L+1] > b
            A[L+1], A[M] = A[M], A[L+1]
            M -= 1
        else # if A[L+1] < b
            L += 1; K += 1
            A[L], A[K] = A[K], A[L]
        end
    end
    return K, M+1
end

#Task_2
function order_statistics!(A, i)
    N = length(A)
    b = A[rand(1:N)]
    K, M = part_sort!(A, b)
    if K < i < M
        return A[i]
    elseif i <= K
        return order_statistics!(@view(A[1:K]), i) 
    else # i >= M
        return order_statistics!(@view(A[M:N]), i)
    end
end
order_statistics(A, i) = order_statistics!(copy(A), i)

#Task_3
function med(A)
    N = length(A)
    if mod(N, 2) == 0
        return (order_statistics!(A, N ÷  2) + order_statistics!(A, N ÷ 2 + 1))/2
    else
        return order_statistics!(A, N/2+1)
    end
end

#Task_4
function minimums(array, k)
    N = length(array)
    k_minimums = sort(array[1:k])
    i = k
    while i < length(array)
        i += 1
        if array[i] < k_minimums[end]
            k_minimums[end] = array[i]
            insert_end!(k_minimums)
        end
    end
    return k_minimums
end         

function insert_end!(array)::Nothing
    j = length(array)
    while j>1 && array[j-1] > array[j]
        array[j-1], array[j] = array[j], array[j-1]
        j -= 1
    end
end

#Task_5
function sd(A)
    s = 0
    s2 = 0
    c = 0
    for i in A
        c+=1
        s+=i
        s2+=i*i
    end
    return sqrt((s2-s^2/c)/c)
end

function comb_sort!(A)
    g = length(A)
    sorted = false
    sm=0
    while !sorted
        g /= 1.237
        if g < 1
            g = 1
            sorted = true
        end
        swapped = false
        for i in 1:length(A)-floor(Int, g)
            sm = floor(Int, g) + i
            if A[i] > A[sm]
                A[i], A[sm] = A[sm], A[i]
                swapped = false
            end
        end
    end
end

f(x) = x

function bubblesort!(a, by::Function = f)
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

A = rand(500)
B = deepcopy(A)

println("bubblesort ", @elapsed bubblesort!(A))
println("comb_sort ", @elapsed comb_sort!(B))
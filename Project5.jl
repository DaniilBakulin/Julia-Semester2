#Task_1
function euclid(m::T, n::T) where T 
    # m, n - заданные целые
    a, b = m, n
    u_a, v_a = T(1), T(0)
    u_b, v_b = T(0), T(1)
    while b != 0
        k = a÷b
        a, b = b, a % b 
        u, v = u_a, v_a
        u_a, v_a = u_b, u_a
        u_b, v_b = u-k*u_b, v-k*v_b
    end
    return a, b, u, v
end

#Task_2 (P7)
function ext_euclid_inv_res(m::T, n::T) where T
    #ax + by = 1    
    a, b = m, n
    u_a, v_a = T(1), T(0)
    u_b, v_b = T(0), T(1)
    a = u_a*m + v_a*n 
    b = u_b*m + v_b*n
    while b != 0
        k = a÷b
        a, b = b, a % b 
        u, v = u_a, v_a
        u_a, v_a = u_b, u_a
        u_b, v_b = u-k*u_b, v-k*v_b
    end
    return u_a
end

#Task_3
function bisect(f::Function, a, b, ε)
    y_a=f(a)
    while b-a > ε
        x_m = (a+b)/2
        y_m=f(x_m)
        if y_m==0
            return x_m
        end
        if y_m*y_a > 0 
            a=x_m
        else
            b=x_m
        end
    end
    return (a+b)/2
end

#Task_4
function isprime_(n::T)::Bool where T <: Integer
    for i in 2:round(T, sqrt(n))
        if n % i == 0
            return false
        end
    end
    return true
end

#Task_5
function eratosphen(n::Integer)::BitVector
    is_prime = ones(Bool, n) 
    is_prime[1] = false # 1 — не является простым числом
    for i in 2:round(Int, sqrt(n))
        if is_prime[i] 
            for j in (i*i):i:n 
                is_prime[j] = false
            end
        end
    end
    return (1:n)[is_prime] # filter(x -> is_prime[x], 1:n) 
end

#Task_6
function dividers(a)
    x = a
    c = 
    d = []
    for i in 2:sqrt(a)
        if i>x
            break
        end
        while mod(x, i) == 0
            c+=1
            x=x/i
        end
        if c!=0
            push!(d, (i, c))
            c=0
        end
    end
    return d
end
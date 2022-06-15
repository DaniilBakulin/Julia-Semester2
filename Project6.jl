#Task_1
function newton(r::Function, x; epsilon = 1e-8, max_num_iter = 20)
    num_iter = 0
    r_x = r(x)
    while num_iter < max_num_iter && abs(r_x) > epsilon
        x += r_x
        r_x = r(x)
        num_iter += 1
    end

    if abs(r_x) <= epsilon
        return x
    else
        return NaN
    end
end

#Task_2
# f(x) = cos(x)-x
# f'(x) = -sin(x)-1
# r(x) = -f(x)/f'(x) = (cos(x)-x)/(sin(x)+1)
newton(x -> (cos(x)-x)/(sin(x)+1), 0.5)

#Task_3
function polyval(P,x)
    dQ = 0
    Q = P[1]
    for i in 2:length(P)        
        dQ = dQ*x + Q
        Q = Q*x + P[i]

    end
    return Q, dQ
end

function r(P, x) #P - внешняя переменная, содержащая коэффициенты многочлена, следующих в порядке убывания степеней
    y, dy = polyval(P, x)
    return -y/dy
end


function complexroot(P::Polynomial, x)
    num_iter = 0
    r_x = r(P, x)
    while num_iter < max_num_iter && abs(r_x) > epsilon
        x += r_x
        r_x = r(P, x)
        num_iter += 1
    end
    if abs(r_x) <= epsilon
        return x
    else
        return NaN
    end
end
#Task_4
function newton2(r::Function, x; epsilon = 1e-8, max_num_iter = 20)
    num_iter = 0
    r_x = -r(x) / ((r(x+h) + r(x)) / epsilon)
    return newton(r=r_x, x)
    
end
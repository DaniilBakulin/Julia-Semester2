using Plots
#import Pkg
#Pkg.add("Plots")

function task_1(n)
    s = 0.0
    a = 1.0
    for k in 1:n+1 
        s += a
        a /= k
    end
    return s
end

function task_1(n)
    s = 0.0
    a = 1.0
    for k in 1:n+1 
        s += a
        a /= k
    end
    return s
end

#Task_2
function eyler(n)
    s = 0.0
    a = 1.0
    for k in 1:n+1 
        s += a
        a /= k
    end
    return s
end

#Task_3
function Base.sin(x,ε)
    xx=x^2
    a=x
    k=1
    s=typeof(x)(0) # - преобразование к 0 нужного типа, что обеспечит стабильность типа переменной s
    while abs(a)>ε
        s+=a
        a=-a*xx/2k/(2k+1)
        k+=1
    end
    #УТВ: |sin(x)-s|<= ε
    return s
end

#Task_4
function exp_(x)
    s = typeof(x)(0) # - преобразование к 0 нужного типа, что обеспечит стабильность типа переменной s
    a = 1
    k = 0
    while s + a != s
        k+=1
        s+=a
        a*=abs(x)
        a/=k
    end
    return s^(sign(x))
end

#Task_5
function harmonic_sum()
    s=0.0
    k=1
    a=1.0
    while s+a != s
        a=1/k
        s+=a 
        k+=1
    end
    return s
end

#Task_6
function cos_n(x, n) #$$\cos(x)=1-\frac{x^2}{2!}+\frac{x^4}{4!}-\frac{x^6}{6!}+...$$
    s = typeof(x)(0)
    a = 1.0
    k=1
    for _ in 1:n+1 
        s += a
        a *= ((-1)*(x^2))
        a /= ((k+1)*k)
        k+=2
    end
    return s
end

#Task_7
function cos_x(x) #$$\cos(x)=1-\frac{x^2}{2!}+\frac{x^4}{4!}-\frac{x^6}{6!}+...$$
    s = typeof(x)(0)
    a = 1.0
    k=1
    while s+a != s
        s += a
        a *= ((-1)*(x^2))
        a /= ((k+1)*k)
        k+=2
    end
    return s
end

#Task_8
function Task_8()
    x = 0:0.1:5
    p = plot()
    for n in 1:5
        plot!(p, x, cos_n.(x, 2^n)) # в текущий график добавлена новая кривая
    end
    display(p) # график, содержащий семейство кривых, отображен
end

#Task_9
function Task_9_1(x) #0 < x < 1
    s = typeof(x)(0)
    a = x-1.0
    k=1
    while s+a != s
        s += a
        a *= ((-1)*(x-1)) 
        a*=k
        k+=1
        a /= k
        
    end
    return s
end

function Task_9_2(x)
    s = typeof(x)(0)
    a = 1
    k=0
    while s+a != s
        s += a
        k+=1
        a *= ((-1)*(2k-1))*x 
        a /= 2k
    end
    return s
end


function Task_9_3(x)
    s = typeof(x)(0)
    a = x^2
    b = x^2
    k=1
    while s+a+b != s
        s += a+b
        k+=1
        a *= -x^2 
        a /= k
        b *= -x^2 
        b /= 2k-2
        b /= 2k-1        
    end
    return s
end

#Task_10
function besselj(m,x)
    s = typeof(x)(0)
    a = 1/factorial(m)
    k=0
    while s+a != s
        s += a
        k += 1
        a *= -(x/2)^2 
        a /= k
        a /= k+m   
    end
    return s*(x/2)^m
end

function Task_10()
    x = 0:0.1:10
    p = plot()
    for m in 0:5
        plot!(p, x, besselj.(m, x)) # в текущий график добавлена новая кривая
    end
    display(p) # график, содержащий семейство кривых, отображен
end
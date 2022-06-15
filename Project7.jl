struct Residue_{T}
    value::T
    M::T
    Residue_{T}(value, M) where T = new(rem(value, M), M) 
end

Base. +(a::Residue_{T}, b::Residue_{T}) where{T} = Residue_{T}(a.value + b.value, a.M)



function Base. *(a::Residue_{T}, b::Residue_{T}) where T
    if a.M==b.M
        return Residue_{T}(a.value * b.value, a.M)
    end
    throw(Exception)
end

Base. -(a::Residue_{T}) where T = Residue_{T}(a.M - a.value, a.M)

function Base. -(a::Residue_{T}, b::Residue_{T}) where T 
    if M.a==M.b
        return Residue_{T}(a.value + a.M-b.value, a.M)
    end
    throw(Exception)
end

function ext_euclid(m::T, n::T) where T
    #ax + by = 1    
    a, b = m, n
    u_a, v_a = 1, 0
    u_b, v_b = 0, 1
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

function Base.inv(a::Residue_{T})::Union{Nothing, Residue_{T}} where T
    if gcd(a.value, M)!=1
        return Nothing
    end

    return Residue_{T}(ext_euclid(a.value, a.M), a.M)
end

function Base. /(a::Residue_{T}, b::Residue_{T}) where T
    if a.M==b.M
        return Residue_{T}(a * inv(b), b.M)
    end
    throw(Exception)
end
Base. ==(a::Residue_{T}, b::Residue_{T}) where T = a.value == b.value
#Base. ===(a::Residue{T,M}, b::Residue{T,M}) where{T,M} = a.value === b.value
Base. >(a::Residue_{T}, b::Residue_{T}) where T = a.value > b.value
Base. <(a::Residue_{T}, b::Residue_{T}) where T = a.value < b.value
Base. <=(a::Residue_{T}, b::Residue_{T}) where T = a.value <= b.value
Base. >=(a::Residue_{T}, b::Residue_{T}) where T = a.value >= b.value

function Ris0(a::Residue_{T}) where T
    return a.value==Int64(0)
end



struct Polynomial{T}
    coeff::Vector{T}
    function Polynomial{T}(coeff) where T 
        n = 0
        if length(coeff) == 1
            return new(coeff)
        end
        for c in reverse(coeff)
            if c != 0
                break
            end
            if c == 0
                n+=1
            end
        end
        new(coeff[1:end-n])
    end
end

deg(p::Polynomial) = length(p.coeff) - 1


function remove_zeros(p::Polynomial{T}) where T
    coeff = copy(p.coeff)
    i, n = lastindex(coeff), 0
    while i > 0 && coeff[i] == 0
        n += 1
        i -= 1
    end
    resize!(coeff, length(coeff)-n)
    return Polynomial{T}(coeff)
end

function Base. +(p::Polynomial{T}, q::Polynomial{T})::Polynomial{T} where T
    np, nq = length(p.coeff), length(q.coeff)
    if  np >= nq 
        coeff = similar(p.coeff)
        coeff[1:nq] .= (@view(p.coeff[1:nq]) .+ q.coeff) 
    else
        coeff = similar(q.coeff)
        coeff[1:np] .= (p.coeff .+ @view(q.coeff[1:np]))
    end
    return remove_zeros(Polynomial{T}(coeff))
end

function Base. -(p::Polynomial{T}, q::Polynomial{T})::Polynomial{T} where T
    np, nq = length(p.coeff), length(q.coeff)
    if  np >= nq 
        coeff = similar(p.coeff)
        coeff[1:nq] .= (@view(p.coeff[1:nq]) .- q.coeff) 
    else
        coeff = similar(q.coeff)
        coeff[1:np] .= (p.coeff .- @view(q.coeff[1:np]))
    end
    # При сложении некоторые старшие коэфициенты могли обратиться в 0 
    return remove_zeros(Polynomial{T}(coeff))
end


function Base. *(p::Polynomial{T}, q::Polynomial{T})::Polynomial{T} where T
    coeff = zeros(T, deg(p) + deg(q)+1)
    for i in eachindex(p.coeff), j in eachindex(q.coeff)
            coeff[i+j - 1] += p.coeff[i]*q.coeff[j]
    end
    return remove_zeros(Polynomial{T}(coeff))
end

Base. +(p::Polynomial{T}, c::T) where T = +(p, Polynomial{T}([c]))
Base. +(c::T, p::Polynomial{T}) where T = +(Polynomial{T}([c]), p)

Base. -(p::Polynomial{T}, c::T) where T = -(p, Polynomial{T}([c]))
Base. -(c::T, p::Polynomial{T}) where T = -(Polynomial{T}([c]), p)

Base. *(p::Polynomial{T}, c::T) where T = *(p, Polynomial{T}([c]))
Base. *(c::T, p::Polynomial{T}) where T = *(Polynomial{T}([c]), p)


function Base.display(p::Polynomial)
    if isempty(p.coeff)
        return ""
    end
    str = "$(p.coeff[1])" # $(...) - означает "интерполяцию стоки", т.е. вставку в строку некоторого вычисляемого значения 
    for i in 2:length(p.coeff)
        if i > 2
            s = " + $(p.coeff[i])x^$(i-1)"
        else
            s = " + $(p.coeff[i])x"
        end
        str *= s
    end
    println(str)
end
(p::Polynomial)(x) = polyval(reverse(p.coeff), x)


function divrem(p::Polynomial{T}, q::Polynomial{T}) where T
    rem = reverse(copy(p.coeff))
    div = []
    for x in 1:deg(p)-deg(q)+1
        println(rem)
        coef = rem[x] / reverse(q.coeff)[1] 
        push!(div, coef)
        for y in 1:deg(q)+1
            rem[x-1+y] = rem[x-1+y] - coef*reverse(q.coeff)[y]
        end
    end
    if length(div) == 0
        div = [0]
    end
    return (Polynomial{T}(div), Polynomial{T}(reverse(rem)))
end
function Base.div(p::Polynomial{T}, q::Polynomial{T})::Polynomial{T} where T 
    return divrem(p,q)[1]
end
function Base.rem(p::Polynomial{T}, q::Polynomial{T})::Polynomial{T} where T 
    return divrem(p,q)[2]
end 


a = Residue_{Int}(5, 2)
c = Polynomial{Int}([1,2,3,4,5,6])
b = Polynomial{Int}([1,2,3])
m = Polynomial{Int}([1,3,3,1])
n = Polynomial{Int}([1,2,1])

A=Residue_{Polynomial{Int}}(c, n)
divrem(m, n)
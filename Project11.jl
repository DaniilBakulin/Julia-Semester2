using LinearAlgebra
using Plots


Vector2D{T<:Real} = NamedTuple{(:x, :y), Tuple{T,T}}

Segment2D{T<:Real} = NamedTuple{(:A, :B), Tuple{Vector2D{T},Vector2D{T}}}

isinner(P::Vector2D, s::Segment2D) = (s.A.x <= P.x <= s.B.x || s.A.x >= P.x >= s.B.x) && 
                                     (s.A.y <= P.y <= s.B.y || s.A.y >= P.y >= s.B.y)


function intersect(s1::Segment2D{T},s2::Segment2D{T}) where T
    A = [s1.B[2]-s1.A[2]  s1.A[1]-s1.B[1]
         s2.B[2]-s2.A[2]  s2.A[1]-s2.B[1]]
    b = [s1.A[2]*(s1.A[1]-s1.B[1]) + s1.A[1]*(s1.B[2]-s1.A[2])
         s2.A[2]*(s2.A[1]-s2.B[1]) + s2.A[1]*(s2.B[2]-s2.A[2])]
    x,y = A\b
    if isinner((;x, y), s1)==false || isinner((;x, y), s2)==false
        return nothing
    end
    return Vector2D{T}((x,y))
end

LinearAlgebra.dot(A::Vector2D{T}, B::Vector2D{T}) where T = dot(Tuple(A), Tuple(B))

function vec_len(A::Vector2D{T}, B::Vector2D{T}) where T
    return sqrt((B[2]-A[2])^2+(B[1]-A[1])^2)
end 


function vec_angle(A::Vector2D{T}, B::Vector2D{T}) where T
    if abs(A.x)==abs(B.x) && abs(A.y)==abs(B.y)
        return acos(dot(A, A)/abs(dot(A,A)))
    end
    return acos(dot(A, B)/(sqrt(dot(A,A)) * sqrt(dot(B,B))))
end 

is_one_area(F::Function, P::Vector2D{T}, Q::Vector2D{T}) where T = (F(P)*F(Q)>0)
function f(X::Vector2D{T}) where T 
    return 2X[1]^2-3X[2]
end

xdot(A::Vector2D{T}, B::Vector2D{T}) where T = A.x*B.y-A.y*B.x

function in_plg(m::Vector2D{T}, plg::Vector{Vector2D{T}}) where T
    v1 = Vector2D{T}((plg[1][1]-m[1], plg[1][2]-m[2]))
    v2 =Vector2D{T}((plg[2][1]-m[1], plg[2][2]-m[2]))
    s = xdot(v1,v2)>0
    for t in 3:length(plg)
        v1 = v2
        v2 = Vector2D{T}((plg[t][1]-m[1], plg[t][2]-m[2]))
        if (xdot(v1, v2)>0)!=s
            return false    
        end
    end
    return true
end

function is_convex(plg::Vector{Vector2D{T}}) where T
    v1 = Vector2D{T}((plg[2][1]-plg[1][1], plg[2][2]-plg[1][2]))
    v2 =Vector2D{T}((plg[3][1]-plg[2][1], plg[3][2]-plg[2][2]))
    s = xdot(v1, v2)>0
    for t in 4:length(plg)
        v1 = v2
        v2 = Vector2D{T}((plg[t][1]-plg[t-1][1], plg[t][2]-plg[t-1][1]))
        if (xdot(v1, v2)>0)!=s
            return false    
        end
    end
    return true
end


function jarvis(plg::Vector{Vector2D{T}}) where T #графики
    startpoint = plg[1]
    wrap = Vector{Vector2D{T}}()
    for i in 1:length(plg)
        if startpoint.y > plg[i].y
            startpoint = plg[i]
        end
    end
    push!(wrap, startpoint)
    nextpoint=0
    nextvector=0
    prevpoint=startpoint
    prevvector = Vector2D{T}((1,0))
    for _ in plg
        minangle = 2pi
        for i in plg
            if i ∉ wrap || i==startpoint#\notin
                vector = Vector2D{T}((i.x-prevpoint.x,i.y-prevpoint.y))
                if minangle > vec_angle(prevvector, vector) && xdot(prevvector, vector)>0
                    nextpoint = i
                    nextvector = vector
                    minangle = vec_angle(prevvector, nextvector)
                end
            end
        end
        if nextpoint==startpoint
            push!(wrap, startpoint)
            return wrap
        end
        prevpoint = nextpoint
        prevvector = nextvector
        push!(wrap, prevpoint)
    end
    return wrap
end

abstract type AbstractPolygon{T<:Real} end

struct Polygon{T} <: AbstractPolygon{T}
    vertices::Vector{Vector2D{T}}
    Polygon{T}(vertices) where T = new(__double_ended!__(vertices))
end

function __double_ended!__(vertices::Vector{Vector2D}) # дублирует в конце вектора его первый элемент, если изначально этого дублирования не было
    if vertices[begin] != vertices[end]
        push!(vertices, polygon[begin])
    end
    return vertices
end

get_vertices(polygon::Polygon) = polygon.vertices
num_vertices(polygon::Polygon) = polygon.vertices[begin] != polygon.vertices[end] ? length(polygon.vertices) : length(polygon.vertices)-1  

struct ConvexPolygon{T} <: AbstractPolygon{T}
    vertices::Vector{Vector2D{T}}
    ConvexPolygon{T}(vertices) where T = new(__double_ended!__(vertices))
end
get_vertices(polygon::ConvexPolygon) = polygon.vertices

function grekhom!(points::Vector{Vector2D{T}})::ConvexPolygon{T} where T<:Real 
    # алгоритм Грехома построения выпуклой оболочки заданного множества точек плоскости
    ydata = (points[i][2] for i in 1:length(points))
    i_start = findmin(ydata) # индекс самой нижней точки
    points[begin], points[i_start] = points[i_start], points[begin]
    sort!(@view(points[begin+1:end]), by = (point -> angle(point, Vector2D{T}(1,0))))
    push!(points, points[begin]) # теперь points[end] == points[begin] 
    convex_polygon = [firstindex(points), firstindex(points)+1, firstindex(points)+2] # - в стек помещены первые 3 точки
    for i in firstindex(points)+3:lastindex(points)
        while sign(points[i]-points[convex_polygon[end]], points[convex_polygon[end-1]]-points[convex_polygon[end]]) < 0
            pop!(convex_polygon)
        end
        push!(convex_polygon, i)
    end
    pop!(points) # из конца массива извлечена предварительно продублированная первая точка (чтобы исходный массив не содержал лишних точек)
    return ConvexPolygon{T}(points[convex_polygon])  # convex_polygon[begin] == convex_polygon[end]
end

function area_tr(plg::Vector{Vector2D{T}}) where T
    prev_p = plg[1]
    s = 0
    for i in 2:length(plg)
        p = plg[i]
        s += (prev_p.y+p.y)*(p.x-prev_p.x)/2
        prev_p = p
    end
    return abs(s)
end


function area_triangle(plg::Vector{Vector2D{T}}) where T
    prev_p = plg[1]
    next_p = plg[1]
    v1 = 0
    v2 = 0
    s = 0
    for i in 2:length(plg)-1
        p = plg[i]
        next_p = plg[i+1]
        v1 = Vector2D{T}((next_p.x - p.x, next_p.y - p.y))
        v2 = Vector2D{T}((prev_p.x - p.x, prev_p.y - p.y))
        s += xdot(v1, v2)
        prev_p = p
    end
    return abs(s)
end

function add_point(plg::Vector{Vector2D{T}}, m::Vector2D{T}) where T
    if in_plg(m, plg)
        return plg
    end
    vector_m = Vector2D{T}((0,0))
    vector_nextp = Vector2D{T}((0, 0))
    vector_currentp = Vector2D{T}((plg[1].x-plg[end-1].x, plg[1].y-plg[end-1].y)) 
    for i in 1:length(plg)-1
        vector_m = Vector2D{T}((m.x-plg[i].x, m.y-plg[i].y))
        vector_nextp = Vector2D{T}((plg[i+1].x-plg[i].x, plg[i+1].y-plg[i].y))
        if (vec_angle(vector_currentp, vector_nextp) >= vec_angle(vector_currentp, vector_m))
            k = i+1
            while true
                if k == length(plg)
                    wrap = plg[:]
                    insert!(wrap, k, m)
                    return wrap
                end
                vector_m = Vector2D{T}((plg[k].x-m.x, plg[k].y-m.y))
                vector_mnext = Vector2D{T}((plg[k+1].x-m.x, plg[k+1].y-m.y))
                println(plg[k], vector_currentp, vector_m, vector_mnext)
                if vec_angle(vector_currentp, vector_m) >= vec_angle(vector_currentp, vector_mnext)
                    println(1)
                    wrap = plg[begin:i]
                    push!(wrap, m) #пересоздается оболчка, затраты памяти
                    append!(wrap, plg[k:end])
                    return wrap
                end
                k+=1
            end
        end
        vector_currentp=vector_nextp
    end
end


function wrap_ap(plg::Vector{Vector2D{T}}) where T
    wrap = plg[1:3]
    push!(wrap, plg[1])
    for p in plg[4:end]
        println(wrap)
        wrap = add_point(wrap, p)
    end
    return wrap
end 

function display_wrap(points::Vector{Vector2D{T}}, wrap::Vector{Vector2D{T}}) where T
    points_x::Vector{T} = []
    points_y::Vector{T} = []
    p = plot()
    for n in wrap
        push!(points_x, n.x)
        push!(points_y, n.y) # в текущий график добавлена новая кривая
    end
    plot!(p, points_x, points_y, legend=false, color="red")
    for n in points
        scatter!(p, [n.x], [n.y], legend=false, color="blue")
    end
    display(p)
end
a=Vector2D{Int}((0,0))
b=Vector2D{Int}((5,0))
c=Vector2D{Int}((5,5))
d=Vector2D{Int}((0,5))
e=Vector2D{Int}((4,5))
o=Vector2D{Int}((5,4))
k=Vector2D{Int}((25,25))
ab = Segment2D{Int}((a, b))
cd = Segment2D{Int}((c, d))
plg1 = Vector{Vector2D{Int}}([b, a, e, k, d, c, o])
plg2 = Vector{Vector2D{Int}}([o, e, d, c, o])
wrap = Vector{Vector2D{Int}}([a, b, o, c, e, d, a])
#add at second place rap_ap([a,o,e,d,c,b])
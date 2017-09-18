
struct Point{T}
    x::T
    y::T
    function Point(x::T,y::T) where T
        @assert length(x) == length(y)
        new{T}(x,y)
    end
end
function show(io::IO,p::Point{T}) where T
    print(io,"Point of type $T:\n")
    print(io,"      x = $(p.x)\n")
    print(io,"      y = $(p.y)\n")
end


mutable struct Line{T<:Number} <: AbstractArray{T<:Number,1}
    x::Vector{T}
    y::Vector{T}
    n::Int
    ex::Tuple
    function Line(x::Vector{T},y::Vector{T}) where {T<:Number}
        this = new{T}()
        this.x = copy(x)
        this.n = length(x)
        this.y = copy(y)
        this.ex = length(x) >0 ? extrema(x) : (0,0)
        @assert length(y)==this.n
        return this
    end
end
function show(io::IO,L::Line{T}) where {T<:Number}
    print(io,"$T Line\n")
    print(io,"number of grid points: $(L.n)\n")
    print(io,"domain: $(L.ex)\n")
    print(io,"range: $(extrema(L.y))\n")
end

function reconfigure!(m::Line)
    # after having updated some objects, need to recompute n
    m.n = length(m.x)
    m.ex = extrema(m.x)
    @assert m.n == length(m.y)
end

eltype(l::Line) = eltype(l.x) 
size(l::Line) = (l.n,)
length(l::Line) = l.n
function getindex(l::Line,i::Int)
    (l.x[i],l.y[i])
end
function getindex(l::Line,i::UnitRange{Int})
    Line(l.x[i],l.y[i])
end
function setindex!(l::Line{T},x::T,y::T,i::Int) where {T<:Number}
    l.x[i] = x
    l.y[i] = y
end
function setindex!(l::Line{T},x::Vector{T},y::Vector{T},i::UnitRange{Int}) where {T<:Number}
    l.x[i] = x
    l.y[i] = y
end

# interpolating a line

function interp(l::Line{T},ix::Vector{T},extrap::Bool=true) where {T<:Number}
    # whenever 
    xex = extrema(ix)
    # @debug(logger,"interpolating $ix ")
    if xex[1] < l.ex[1] || xex[2] > l.ex[2]
        if extrap 
            itp = extrapolate(interpolate((l.x,),l.y,Gridded(Linear())),Linear())
        else
            itp = extrapolate(interpolate((l.x,),l.y,Gridded(Linear())),-Inf)
        end
    else
        itp = interpolate((l.x,),l.y,Gridded(Linear()))
    end
    return itp[ix]
end 

# appending, prepending , deleting and splitting at

"prepend points to an Mline"
function prepend!(m::Line,x,y)
    prepend!(m.x,x)
    prepend!(m.y,y)
    reconfigure!(m)
end

"delete an index"
function delete!(m::Line,idx::Int)
    deleteat!(m.x,idx)
    deleteat!(m.y,idx)
    reconfigure!(m)
end

"append points to an Mline"
function append!(m::Line,x,y)
    append!(m.x,x)
    append!(m.y,y)
    reconfigure!(m)
end

"insert a single value at an interior index"
function insert!(m::Line,vx,vy,idx::Int ) 
    insert!(m.x,idx,vx)
    insert!(m.y,idx,vy)
    reconfigure!(m)
end

"""
    splitat(m::Line,j::Int,repeat_boundary::Bool=true)

Splits a `Line` object after given index and returns 2 new `Line`s as a tuple. If `repeat_boundary` is true, then the separating index is the first point of the second new `Line`.
"""
function splitat(m::Line,j::Int,repeat_boundary::Bool=true)
    m1 = Line(m.x[1:j],m.y[1:j])
    if repeat_boundary
        m2 = Line(m.x[j:end],m.y[j:end])
    else
        m2 = Line(m.x[j+1:end],m.y[j+1:end])
    end
    return (m1,m2)
end

"sort a `Line` along x-grid"
function sort!(m::Line)
    ix = sortperm(m.x)
    m.x = m.x[ix]
    m.y = m.y[ix]
end


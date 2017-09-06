
module MultiLine
using StaticArrays
using Interpolations

export Line, interp, append!, prepend!, insert!, delete!

import Base.size, 
       Base.getindex, 
       Base.setindex!, 
       Base.eltype


# started with this, because ranges are super efficient with interpolations, but too inflexible for this purpose.
# mutable struct Line{T<:Number} <: AbstractArray{T<:Number,1}
#     x::Union{Vector{T},StepRangeLen{T}}
#     y::Vector{T}
#     n::Int
#     function Line(x::Union{Vector{T},StepRangeLen{T}},y::Vector{T}) where {T<:Number}
#         this = new{T}()
#         this.x = x
#         this.n = length(x)
#         this.y = y
#         @assert length(y)==this.n
#         return this
#     end
# end
mutable struct Line{T<:Number} <: AbstractArray{T<:Number,1}
    x::Vector{T}
    y::Vector{T}
    n::Int
    function Line(x::Vector{T},y::Vector{T}) where {T<:Number}
        this = new{T}()
        this.x = x
        this.n = length(x)
        this.y = y
        @assert length(y)==this.n
        return this
    end
end
function reconfigure!(m::Line)
    # after having updated some objects, need to recompute n
    m.n = length(m.x)
    @assert m.n = length(m.y)
end

eltype(l::Line) = eltype(l.x) 
size(l::Line) = (l.n,)
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

function interp(l::Line{T},ix::T) where {T<:Number}
    itp = Interpolations.interpolate((l.x,),l.y,Gridded(Linear()))
    return itp[ix]
end 
# function interp(x::StepRangeLen{T},y::Vector{T},ix::T) where {T<:Number}

#     itp = Interpolations.interpolate(y,Bspline(Linear()))
#     sitp = scale(itp,x)
#     return sitp[ix]
# end

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


"""
    secondary_envenlope(m::Mline)

Prunes the `Mline` object from wrong EGM solution points. Wrong solutions appear in kinked regions.
"""
function secondary_envelope(m::Line)
    o = copy(m)  # make a copy in order not to destroy the original object m. check if it's necessary to keep that object, otherwise: destroy it and don't copy anything!

    # 1) find all jump-backs in x-grid
    ii = o.x[2:end].>o.x[1:end-1]  

    # 2) if no backjumps at all, exit
    if all(ii)  
        return m
    else
    # 3) else, identify subsets
        i = 1
        sections = Line[]  # an array of Mlines
        while true
            j = findfirst(ii .!= ii[1])  # identifies all indices within kinked region from left to right until the first kink

            # if no more kinks
            if j==0
                if i > 1
                    # add remaining Mline
                    push!(sections,o)
                end
                break
            end
            newm,o = splitat(o,j)  # split old Mline at j 
            push!(sections,newm)
            ii = ii[j:end] # chop off from jump index
            i += 1
        end

        # 4) do secondary envelope pruning
        # sort all sections on x
        for s in sections
            sort!(s)
        end
        # 5) compute upper envelope of all sections
            # - get all x's from all s and sort into a vector xx
            # - interpolate(extrapolate) all s on xx
            # - disregard points where some lines are extrapolated


        # 6) collect indices of points removed from o (i.e. points in m but not in o)
    end
end

function upper_env(m::Line)
    # 5) compute upper envelope of all sections
        # - get all x's from all s and sort into a vector xx
        # - interpolate(extrapolate) all s on xx
        # - disregard points where some lines are extrapolated


end



end # module



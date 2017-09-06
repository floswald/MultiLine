
module MultiLine
using StaticArrays
using Interpolations

export Mline, interpolate

"""
    multiple interpolations object

Stores a one-dimensional function support ``x`` as well as corresponding
function values ``y``, where ``y`` can have more than one dimensions.

## Fields

* `x`: support of function on 1-dimensional grid
* `y`: function values, `m`-dimensional
* `m`: number of dimensions of `y`
* `n`: number of points in `x` and `y` grids
* `mn`: total number of points, i.e. `n` times `m`

"""
mutable struct Mline 
    x :: Vector{Float64}
    y :: Array{Vector{Float64},1}
    m :: Int
    n :: Int
    nm ::Int
    function Mline(x::Vector{Float64},y::Array{Vector{Float64},1})
        this = new()
        this.m = length(y)
        this.n = length(x)
        @assert all( [this.n == length(i) for i in y] )
        this.x = x
        this.y = y
        this.nm = this.n*this.m
        return this
    end
end
function reconfigure!(m::Mline)
    # after having updated some objects, need to recompute n and am
    m.n = length(m.x)
    m.m = length(m.y)
    m.nm = m.n*m.m
end

function copy(m::Mline)
    mm = Mline(deepcopy(m.x),deepcopy(m.y))
    return mm
end

# methods that change an Mline object

"delete an index"
function delete!(m::Mline,idx::Int)
    deleteat!(m.x,idx)
    for i in m.y
        deleteat!(i,idx)
    end
    reconfigure!(m)
end

"append points to an Mline"
function append!(m::Mline,x,y::Vector)
    append!(m.x,x)
    for i in eachindex(m.y)
        append!(m.y[i],y[i])
    end
    reconfigure!(m)
end

"prepend points to an Mline"
function prepend!(m::Mline,x,y::Vector)
    prepend!(m.x,x)
    for i in eachindex(m.y)
        prepend!(m.y[i],y[i])
    end
    reconfigure!(m)
end

"insert a single value at an interior index"
function insert!(m::Mline,vx::Float64,vy::Vector{Float64},idx::Int ) 

    insert!(m.x,idx,vx)
    for i in eachindex(m.y)
        insert!(m.y[i],idx,vy[i])
    end
    reconfigure!(m)
end

"insert multiple values at an index"
function insert!(m::Mline,vx::Vector{Float64},vy::Vector{Vector{Float64}},idx::Int=1 ) 
    # this needs to create new vectors
    m.x = vcat()
    for i in eachindex(m.y)
        insert!(m.y[i],idx,vy[i])
    end
    reconfigure!(m)
end

"sort an `Mline` along x-grid"
function sort!(m::Mline)
    ix = sortperm(m.x)
    m.x = m.x[ix]
    for i in eachindex(m.y)
        m.y[i] = m.y[i][ix]
    end
end

"""
    splitat(m::Mline,j::Int,repeat_boundary::Bool=true)

Splits an `Mline` object after given index and returns 2 new `Mline`s as a tuple. If `repeat_boundary` is true, then the separating index is the first point of the second new `Mline`.
"""
function splitat(m::Mline,j::Int,repeat_boundary::Bool=true)
    m1 = Mline(m.x[1:j],m.y[1:j])
    if repeat_boundary
        m2 = Mline(m.x[j:end],m.y[j:end])
    else
        m2 = Mline(m.x[j+1:end],m.y[j+1:end])
    end
    return (m1,m2)
end

# methods that produce a value from an Mline object.

"interpolate at new index ix"
function interpolate(m::Mline,ix::Float64)
    if m.m > 1
        # cast as a m-dimensional static array
        y = reinterpret(SVector{m.m,Float64},vcat(m.y'...),(m.n,))
        itp = Interpolations.interpolate(y,BSpline(Linear()), OnGrid())
    else
        itp = Interpolations.interpolate(m.y,BSpline(Linear()), OnGrid())
    end
    return itp[ix]
end

"""
    secondary_envenlope(m::Mline)

Prunes the `Mline` object from wrong EGM solution points. Wrong solutions appear in kinked regions.
"""
function secondary_envelope(m::Mline)
    o = copy(m)  # make a copy in order not to destroy the original object m. check if it's necessary to keep that object, otherwise: destroy it and don't copy anything!

    # 1) find all jump-backs in x-grid
    ii = o.x[2:end].>o.x[1:end-1]  

    # 2) if no backjumps at all, exit
    if all(ii)  
        return m
    else
    # 3) else, identify subsets
        i = 1
        sections = Mline[]  # an array of Mlines
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


function test()
    x = collect(0:0.1:1)
    y = Vector{Float64}[]
    push!(y,rand(11),rand(11),rand(11))
    m = Mline(x,y)
    interpolate(m,3.1234)
end


end # module



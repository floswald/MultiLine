
module MultiLine

using StaticArrays
using Interpolations
using Roots
using MiniLogging

export Line, Point, interp, splitat,upper_env

import Base.size, 
       Base.getindex, 
       Base.setindex!, 
       Base.eltype,
       Base.prepend!,
       Base.append!,
       Base.insert!,
       Base.delete!,
       Base.sort!

# setup MiniLogging
logger = get_logger()
# basic_config(MiniLogging.DEBUG; date_format="%Y-%m-%d %H:%M:%S")
basic_config(MiniLogging.INFO; date_format="%Y-%m-%d %H:%M:%S")


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

struct Point{T}
    x::T
    y::T
end

mutable struct Line{T<:Number} <: AbstractArray{T<:Number,1}
    # pts::Vector{Point{T}}
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
function Line() 
    Line(Number[],Number[])
end
function reconfigure!(m::Line)
    # after having updated some objects, need to recompute n
    m.n = length(m.x)
    m.ex = extrema(m.x)
    @assert m.n == length(m.y)
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

function interp(l::Line{T},ix::Vector{T},extrap::Bool=true) where {T<:Number}
    # whenever 
    xex = extrema(ix)
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

"Interpolate a Vector of `Line`s on the same grid. Return a matrix where each row is the interpolation of another `Line`"
function interp(L::Vector{Line{T}},ix::Vector{T},extrap::Bool=true) where {T<:Number}

    # yy = reinterpret(SVector{length(L),T},vcat([l.y for l in L]'...),(L[1].n,))
    # itp = interpolate((xx,),yy,Gridded(Linear()))
    # y = itp[ix]
    # y = convert(Matrix{T},y)
    y = zeros(T,length(L),length(ix))
    for i in eachindex(L)
        y[i,:] = interp(L[i],ix,extrap)
    end
    return y
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
function secondary_envelope(o::Line{T}) where T<:Number

    # 1) find all jump-backs in x-grid
    ii = o.x[2:end].>o.x[1:end-1]  

    # 2) if no backjumps at all, exit
    if all(ii)  
        # return same object
        return Dict(:pruned=>o,:idx_removed=>0,:intersections=>0)
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
                # then break
                break
            end
            newm,o = splitat(o,j)  # split old Line at j 
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
        u_env = upper_env(s)
        return Dict(:pruned=>u_env[:envelope],:idx_removed=>diff(o,u_env[:envelope]),:intersections=>u_env[:intersections])
    end
end

function upper_env(L::Vector{Line{T}}) where T<:Number
    # 5) compute upper envelope of all lines
        # - get all x's from all s and sort into a vector xx
        # - interpolate(extrapolate) all s on xx
        # - how to deal with points at which some Line is infeasible?

    if length(L)==1
        warn("an upper envelope requires by definition at least 2 lines.")
        return 0
    end

    # - get all x's from all Lines and sort into a vector xx
    xx = sort(unique(vcat([l.x for l in L]...)))
    n = length(xx)

    # - interpolate(extrapolate) all Ls on xx
    # this returns a matrix (length(L),n)
    # i.e. each row is the interpolation 
    yy = interp(L,xx)

    # find the top line at each point in xx
    val,lin_ind = findmax(yy,1)  # colwise max
    subs = map(x->ind2sub(yy,x),lin_ind)  # get subsript indices of colwise maxima
    r_idx = [i[1] for i in subs] # get row indices only: the row index tells us which Line was optimal at that point.

    # Identify changes in optimal Line
    # switch in top line after index s (indexing global support xx)
    # s tells us after which position in xx we have a change in optimal line
    s = find(r_idx[2:end].!=r_idx[1:end-1])


    # Assemble Upper Envelope from Line segments
    # ==========================================

    if length(s)==0
        # there is one complete upper envelope already
        # return
        return Dict(:envelope => yy[r_idx[1],:],:intersections=>0)
    else
        # sort out which line is top at which index of xx and compute intersection points in between switches
        # s = 1: there is a switch in top line after the first index
        # s = i: there is a switch in top line after the i-th index

        # return Line 
        # The envelope starts with the first line that is on top
        # that line is on top until index s[1] in xx, after which the 
        # top line changes.
        env = Line(xx[1:s[1]],yy[r_idx[s[1]],1:s[1]] )
        isec = Point{T}[]

        for id_s in eachindex(s)

            js = s[id_s]  # value of index: position in xx
            # @debug(logger,"js = $js")

            # switching from Line to Line
            from = r_idx[js]
            to   = r_idx[js+1]
            @debug(logger,"from = $(r_idx[js])")
            @debug(logger,"to   = $(r_idx[js+1])")

            # xx coordinates between which the switching occurs
            # remember xx is a vector as long as size(yy,2)
            x_from = xx[subs[js][2]]  # only pick col coordinate
            x_to   = xx[subs[js+1][2]]  # only pick col coordinate
            @debug(logger,"x_from = $(xx[subs[js][2]])")  # only pick col coordinate
            @debug(logger,"x_to   = $(xx[subs[js+1][2]])")  # only pick col coordinate

            # end and start values of both lines
            v_from = yy[subs[js]...]
            v_to   = yy[subs[js+1]...]
            @debug(logger,"v_from = $(yy[subs[js]...])")
            @debug(logger,"v_to   = $(yy[subs[js+1]...])")

            # if both L[from] and L[to] have an xgrid such that
            # xx=[...,x_from,x_to,...] both x_from and x_to are members,
            # then there is no intersection to compute: the switch happens on a 
            # grid point, hence the intersection is x_to. 
            # also, we don't have to add the intersection to the envelope
            # (because x_to would show up twice in this case: last point of L[from] and
            # intersection).
            from_on_xx = to_on_xx = false

            # check if L[from].x contains adjacent x_from and x_to
                # @debug(logger,in(x_from,L[from].x))
                # @debug(logger,L[from].x[findfirst(L[from].x,x_from)+1])
                # @debug(logger,x_to)
            if in(x_from,L[from].x) && L[from].x[findfirst(L[from].x,x_from)+1]==x_to
                from_on_xx = true
            end
            if in(x_from,L[to].x) && L[to].x[findfirst(L[to].x,x_from)+1]==x_to
                to_on_xx = true
            end
            @debug(logger,"from_on_xx = $from_on_xx")
            @debug(logger,"to_on_xx = $to_on_xx")

            # if both lines' support is contained in global support xx
            if from_on_xx && to_on_xx
                # record intersection
                push!(isec,Point(x_to,v_to))

                # don't add intersection to envelope!

                # add next line segment to envelope
                # index range s[id_s]+1:s[id_s+1] is
                #     from current switch (next index after): s[id_s]+1
                #     to last index before next switch: s[id_s+1]
                last_ind = id_s==length(s) ? n : s[id_s+1]
                append!(env,xx[js+1:last_ind],yy[to,js+1:last_ind])
            else
                # compute location in grid and value at intersection
                # @debug(logger,interp(L[to],[x_to-x_from]))
                # @debug(logger,interp(L[from],[x_to-x_from]))
                # @debug(logger,x_from)
                # @debug(logger,x_to)
                f_closure(x) = interp(L[to],[x])[1] - interp(L[from],[x])[1]
                x_x = fzero(f_closure, x_from, x_to)
                v_x = interp(L[to],[x_x])[1]

                # record intersection
                push!(isec,Point(x_x,v_x))
                @debug(logger,"isec = $isec")

                # if intersection is not a grid point already,
                # add intersection to envelope
                # @debug(logger,"in(x_x,xx) = $(in(x_x,xx))")
                # @debug(logger,"maxdiff(x_x,xx) = $(minimum(abs,(x_x-xx)))")
                # @debug(logger,"maxdiff(v_x,xx) = $(minimum(abs,(v_x-yy)))")
                if !(minimum(abs.(x_x-xx)) < sqrt(eps()) && minimum(abs.(v_x-yy)) < sqrt(eps())  )
                    append!(env,x_x,v_x)
                end

                # add next line segment to envelope
                # index range s[id_s]+1:s[id_s+1] is
                #     from current switch (next index after): s[id_s]+1
                #     to last index before next switch: s[id_s+1]
                last_ind = id_s==length(s) ? n : s[id_s+1]
                append!(env,xx[js+1:last_ind],yy[to,js+1:last_ind])
            end

        end
        return Dict(:envelope => env,:intersections=>isec)
    end
end



end # module



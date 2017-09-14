
"""

## Envelope

Holds an array of `Line`s, the upper envelope of those lines, and a vector of `Point`s marking the intersections between lines.
"""
mutable struct Envelope{T<:Number}
    L      :: Vector{Line{T}}
    env    :: Line{T}
    isects :: Vector{Point{T}}
    removed :: Vector{Vector{Point{T}}}
    # function Envelope(l::Vector{Line{T}},e::Vector{T},ise::Vector{Point{T}}) where {T<:Number}
    #     this = new{T}()
    #     this.L = l
    #     this.env = e 
    #     this.isects = ise
    #     this.removed = [Point{T}[]]
    #     return this
    # end
    function Envelope(e::Line{T}) where {T<:Number}
        this = new{T}()
        this.L = Line{T}[]
        this.env = e 
        this.isects = Point{T}[]
        this.removed = Vector{Point{T}}[Point{T}[] ]
        return this
    end
    function Envelope(l::Vector{Line{T}}) where {T<:Number}
        this = new{T}()
        this.L = l
        this.env = Line([typemin(T)],[typemin(T)])
        this.isects = Point{T}[]
        this.removed = Vector{Point{T}}[Point{T}[] ]
        return this
    end
end
function show(io::IO,en::Envelope{T}) where {T<:Number}
    print(io,"$T Envelope\n")
    print(io,"num of `Line`s: $(length(en.L))\n")
    print(io,"num of intersections: $(length(en.isects))\n")
    print(io,"num of pts removed: $(length(en.removed))\n")
end
size(e::Envelope) = size(e.L)
eltype(e::Envelope) = eltype(e.L) 
getx(en::Envelope) = en.env.x
gety(en::Envelope) = en.env.y
gets(en::Envelope) = en.isects
getr(en::Envelope) = en.removed



"Interpolate an `Envelope` on a unique grid. Return a matrix where each row is the interpolation of another `Line`"
function interp(e::Envelope{T},ix::Vector{T},extrap::Bool=true) where {T<:Number}

    # yy = reinterpret(SVector{length(L),T},vcat([l.y for l in L]'...),(L[1].n,))
    # itp = interpolate((xx,),yy,Gridded(Linear()))
    # y = itp[ix]
    # y = convert(Matrix{T},y)
    y = zeros(T,length(e.L),length(ix))
    for i in eachindex(e.L)
        y[i,:] = interp(e.L[i],ix,extrap)
    end
    return y
end 

# function interp(x::StepRangeLen{T},y::Vector{T},ix::T) where {T<:Number}

#     itp = Interpolations.interpolate(y,Bspline(Linear()))
#     sitp = scale(itp,x)
#     return sitp[ix]
# end

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



"""
    create_envelope(m::Mline)

splits a `Line` object at wrong EGM solution points. Wrong solutions appear in kinked regions.
"""
function create_envelope(o::Line{T}) where T<:Number

    # 1) find all jump-backs in x-grid
    ii = o.x[2:end].>o.x[1:end-1]  

    # 2) if no backjumps at all, exit
    if all(ii)  
        # return as an Envelope
        return Envelope(o)
    else
    # 3) else, identify subsets
        i = 1
        sections = Line{T}[]  # an array of Lines
        while true
            # println(ii)
            j = findfirst(ii .!= ii[1])  # identifies all indices within kinked region from left to right until the first kink
            # println(j)

            # if no more kinks
            if j==0
                if i > 1
                    # add remaining Line
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

        # 4) sort all sections on x
        for s in sections
            sort!(s)
        end
        return Envelope(sections)
    end
end

function upper_env!(e::Envelope{T}) where T<:Number
    # 5) compute upper envelope of all lines
        # - get all x's from all s and sort into a vector xx
        # - interpolate(extrapolate) all s on xx
        # - how to deal with points at which some Line is infeasible?

    if length(e.L)<2
        error("an upper envelope requires by definition at least 2 lines.")
    end

    # - get all x's from all Lines and sort into a vector xx
    xx = sort(unique(vcat([l.x for l in e.L]...)))
    n = length(xx)

    # - interpolate(extrapolate) all Ls on xx
    # this returns a matrix (length(L),n)
    # i.e. each row is the interpolation 
    yy = interp(e,xx)

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
        e.env = Line(xx,yy[r_idx[1],:])
        e.isects = [Point(NaN,NaN)]
        e.removed = [[Point(NaN,NaN)]]
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
            jx  = subs[js][2]   # colindex of element in yy before switch takes place
            jjx = subs[js+1][2] # colindex of element in yy after switch took place

            # switching from Line to Line
            from = r_idx[js]
            to   = r_idx[js+1]
            @debug(logger,"from = $(r_idx[js])")
            @debug(logger,"to   = $(r_idx[js+1])")

            # xx coordinates between which the switching occurs
            # remember xx is a vector as long as size(yy,2)
            x_from = xx[jx ]  # only pick col coordinate
            x_to   = xx[jjx]  # only pick col coordinate
            @debug(logger,"x_from = $(xx[jx ])")  # only pick col coordinate
            @debug(logger,"x_to   = $(xx[jjx])")  # only pick col coordinate

            # values of both lines at those corrdinates
            v_from = yy[subs[js]...]
            v_to   = yy[subs[js+1]...]
            @debug(logger,"v_from = $(yy[subs[js]...])")
            @debug(logger,"v_to   = $(yy[subs[js+1]...])")

            # The intersction of both lines ∈ [x_from,x_to]
            # complication: intersection could also be on the egdes of this interval.
            # this is easy to check:
            if isapprox(yy[from,jx] , yy[to,jx])
                # intersection is on lower bound of interval
                push!(isec,Point(x_from,v_from))
                # don't add intersection to envelope!
            elseif isapprox(yy[from,jjx] , yy[to,jjx])
                # intersection is on upper bound of interval
                push!(isec,Point(x_to,v_to))
                # don't add intersection to envelope!
            else
                # need to to compute intersection
                f_closure(z) = interp(e.L[to],[z])[1] - interp(e.L[from],[z])[1]
                x_x = fzero(f_closure, x_from, x_to)
                v_x = interp(e.L[from],[x_x])[1]

                # record intersection
                push!(isec,Point(x_x,v_x))
                # and add to envelope
                append!(env,x_x,v_x)
                @debug(logger,"added $(Point(x_x,v_x)) to envelope")
                @debug(logger,"x in env? $(in(x_x,env.x))")
            end

            # add next line segment to envelope
            # index range s[id_s]+1:s[id_s+1] is
            #     from current switch (next index after): s[id_s]+1
            #     to last index before next switch: s[id_s+1]
            last_ind = id_s==length(s) ? n : s[id_s+1]
            append!(env,xx[js+1:last_ind],yy[to,js+1:last_ind])
        end
        e.env = env 
        e.isects = isec
        # collect points that were removed from Lines
        for l in e.L
            ix = find( !in(getx(e),l.x) || !in(gety(e),l.y) )
            push!(e.removed,[Point(l.x[jx],l.y[jx]) for jx in ix])
        end
        return nothing
    end
end



module MultiLine

# using StaticArrays
using Interpolations
using Roots
using Plots
# plotlyjs()
pyplot()

# Types
export Line, Point, Envelope

# methods
export interp, splitat,upper_env!, getx, gety, gets, create_envelope, getr

import Base.size, 
       Base.getindex, 
       Base.setindex!, 
       Base.eltype,
       Base.prepend!,
       Base.append!,
       Base.insert!,
       Base.delete!,
       Base.sort!,
       Base.length,
       Base.show


# includes
include("line.jl")
include("envelope.jl")
include("plotting.jl")


end # module



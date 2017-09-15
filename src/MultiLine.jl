
module MultiLine

using StaticArrays
using Interpolations
using Roots
using MiniLogging

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


# setup MiniLogging
logger = get_logger()
# if isinteractive()
    basic_config(MiniLogging.DEBUG; date_format="%H:%M:%S")
# else
#     basic_config(MiniLogging.INFO; date_format="%H:%M:%S")
# end


# includes
include("line.jl")
include("envelope.jl")


end # module



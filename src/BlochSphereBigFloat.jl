module BlochSphereBigFloat

using StaticArrays
using CairoMakie
using ProgressMeter
using CoordinateSystems
import Base: show, convert, propertynames, getproperty, iterate, length

export SingleQubitState, 
    BlochVec,

    StereographicCoordinates,
    projection,

    setup_blochplot,
    setup_doubleplot,
    setup_tripleplot,
    plot_axis_arrow,
    interpolate_path,

    PauliX, PauliY, PauliZ,
    ep, Rx, Ry, Rz

include("SingleQubitState.jl")
include("stereographic-projection.jl")
include("visualization.jl")
include("helpers.jl")

include("pulse-sequence-randomized-rotations.jl")
#include("video-fix.jl")




end # module BlochSphereBigFloat

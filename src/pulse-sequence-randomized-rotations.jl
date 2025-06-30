# submodule of BlochSphereBigFloat
"""
# submodule PulseSequenceRandomizedRotations

Contains functions to animate and calculate states using pulse sequences with arbitrary x- and randomized z-rotations.
"""
module PulseSequenceRandomizedRotations

using ThreadsX
using Base.Threads
using StaticArrays
using ProgressMeter
using ComplexBigMatrices
using LinearAlgebra
using CairoMakie
using CoordinateSystems
using BlochSphereBigFloat

export animate_sequence_path, 
    animate_sequence_path_triple,
    plot_sequence_path,
    plot_sequence_path_triple,

    states_by_sequence,
    overlap_by_sequence,
    states_by_sequence_threaded,
    overlap_by_sequence_threaded,
    vals2sequence,
    seq2vals


# Setup functions for plotting
function plot_sequence_path(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}, N::Int64; s_start=SingleQubitState([1,1]), pBlochargs=(colormap=:viridis,), kwargs...)
    fig, ax_bloch = setup_blochplot(;kwargs...)
    dt = 1/N
    
    states = Vector[SingleQubitState[s_start] for _ in ωs]
    n = length(ωs)
    
    for pulse in pulses
        dir, t = pulse
        for i in 1:n
            if dir == :z
                Rot = Rz
                ϕ = ωs[i] * t
            else
                Rot = Rx
                ϕ = t
            end
            push!(states[i], interpolate_path(states[i][end], ϕ, Rot, N)...)
        end
    end
    
    for i in 1:n
        pBloch.(ax_bloch, states[i], false, color=i, colormap=:viridis, colorrange=(1,n+1), alpha=0.7; pBlochargs...)
        pBloch(ax_bloch, states[i][end], true, color=i, colormap=:viridis, colorrange=(1,n+1); pBlochargs...)
    end
    
    fig, ax_bloch
    
end

function plot_sequence_path_triple(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}, N::Int64; kwargs...)
    fig, ax_bloch, ax_stereo, ax_stereo_s = setup_tripleplot(;kwargs...)
    dt = 1/N
    
    s_plus = SingleQubitState([1,1])
    
    states = Vector[SingleQubitState[s_plus] for _ in ωs]
    n = length(ωs)
    
    for pulse in pulses
        dir, t = pulse
        for i in 1:n
            if dir == :z
                Rot = Rz
                ϕ = ωs[i] * t
            else
                Rot = Rx
                ϕ = t
            end
            push!(states[i], interpolate_path(states[i][end], ϕ, Rot, N)...)
        end
    end
    
    for i in 1:n
        pBloch.(ax_bloch, states[i], false, color=i, colormap=:viridis, colorrange=(1,n+1), alpha=0.7)
        pBloch(ax_bloch, states[i][end], true, color=i, colormap=:viridis, colorrange=(1,n+1))
        states_stereo = stack(convert.(StereographicCoordinates{PolarCoordinates}, projection.(StereographicCoordinates, states[i])))
        scatter!(ax_stereo, states_stereo[2,:], states_stereo[1,:], color=i, colormap=:viridis, colorrange=(1,n+1), alpha=0.7) 
        states_stereo_s = stack(convert.(StereographicCoordinates{PolarCoordinates}, projection.(StereographicCoordinates, states[i], proj=:south)))
        scatter!(ax_stereo_s, states_stereo_s[2,:], states_stereo_s[1,:], color=i, colormap=:viridis, colorrange=(1,n+1), alpha=0.7) 
    end
    
    fig, ax_bloch, ax_stereo, ax_stereo_s
    
end

function animate_sequence_path(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}, N::Int64; s_start=false, plot_trace=true, plot_step_trace=false, plot_arrow=false, fps=24, markersize=5, kwargs...)
    fig, ax_bloch = setup_blochplot(;kwargs...)
    dt = 1/N
    
    if s_start == false
        s_start = SingleQubitState([1,1]) # |+>
    end
    
    states = Vector[SingleQubitState[s_start] for _ in ωs]
    n = length(ωs)
    
    # generate all states
    for pulse in pulses
        dir, t = pulse
        for i in 1:n
            if dir == :z
                Rot = Rz
                ϕ = ωs[i] * t
            else
                Rot = Rx
                ϕ = t
            end
            push!(states[i], interpolate_path(states[i][end], ϕ, Rot, N)...)
        end
    end
    
    # overlaps
    s_z = SingleQubitState([1,0])
    overlaps = [abs2(dot(state[end], s_z)) for state in states]
    cmap = haskey(kwargs, :colormap) ? kwargs[:colormap] : :viridis

    plotobjs_arrow = []
    plotobjs_trace = []

    # animate
    anim = Record(fig, framerate=fps) do io
        @showprogress "Animating..." for j in 1:length(states[1])
            if !plot_trace && !plot_step_trace
                delete!.(ax_bloch, plotobjs_trace)
            elseif plot_step_trace && length(plotobjs_trace) >= n*N
                delete!.(ax_bloch, plotobjs_trace)
                plotobjs_trace = [] # clear all traces
            end

            if plot_arrow
                delete!.(ax_bloch, plotobjs_arrow)
            end
            for i in 1:n
                if !plot_trace && !plot_step_trace # no trace
                    p = pBloch(ax_bloch, states[i][j], plot_arrow, color=i, colormap=cmap, colorrange=(1,n+1), alpha=0.7)
                    push!(plotobjs_trace, p)
                elseif plot_arrow # trace and arrow, which is the last state
                    # trace
                    p1 = pBloch(ax_bloch, states[i][j], false, color=i, colormap=cmap, colorrange=(1,n+1), alpha=0.7, markersize=markersize)
                    # arrow
                    p2 = pBloch(ax_bloch, states[i][j], true, color=i, colormap=cmap, colorrange=(1,n+1), alpha=0.7)
                    # store the arrow for deletion in the next frame
                    push!(plotobjs_trace, p1)
                    push!(plotobjs_arrow, p2)
                else # only trace
                    pBloch(ax_bloch, states[i][j], false, color=i, colormap=cmap, colorrange=(1,n+1), alpha=0.7, markersize=markersize)
                end
            end
            recordframe!(io)
        end
    end
    
    anim
    
end

function animate_sequence_path_triple(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}, N::Int64; s_start=false, plot_trace=true, plot_arrow=false, rmax=1.5, fps=24, kwargs...)
    fig, ax_bloch, ax_stereo, ax_stereo_s = setup_tripleplot(;kwargs...)
    rlims!(ax_stereo, 0,rmax)
    rlims!(ax_stereo_s, 0,rmax)
    dt = 1/N
    
    if s_start == false
        s_start = SingleQubitState([1,1]) # |+>
    end
    
    states = Vector[SingleQubitState[s_start] for _ in ωs]
    n = length(ωs)
    
    # generate all states
    for pulse in pulses
        dir, t = pulse
        for i in 1:n
            if dir == :z
                Rot = Rz
                ϕ = ωs[i] * t
            else
                Rot = Rx
                ϕ = t
            end
            push!(states[i], interpolate_path(states[i][end], ϕ, Rot, N)...)
        end
    end

    plotobjs = []

    # animate
    anim = Record(fig, framerate=fps) do io
        @showprogress "Animating..." for j in 1:length(states[1])
            if !plot_trace
                delete!.(ax_bloch, plotobjs)
            end
            for i in 1:n
                p = pBloch(ax_bloch, states[i][j], !plot_trace && plot_arrow, color=i, colormap=:viridis, colorrange=(1,n+1), alpha=0.7)
                stereo = convert(StereographicCoordinates{PolarCoordinates}, projection(StereographicCoordinates, states[i][j])) 
                scatter!(ax_stereo, [stereo.azi], [stereo.r], color=i, colormap=:viridis, colorrange=(1,n+1), alpha=0.7) 
                stereo_s = convert(StereographicCoordinates{PolarCoordinates}, projection(StereographicCoordinates, states[i][j], proj=:south)) 
                scatter!(ax_stereo_s, [stereo_s.azi], [stereo_s.r], color=i, colormap=:viridis, colorrange=(1,n+1), alpha=0.7) 
                if !plot_trace
                    push!(plotobjs, p)
                end
            end
            recordframe!(io)
        end
    end
    
    anim
    
end

# Basic functions for pulse sequence randomized rotations
# BASIC NEEDED FUNCTIONS
"""
    function states_by_sequence(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1))

Calculates the states of a single qubit after applying a sequence of pulses defined by `pulses` to the initial state `s_start`. Each pulse is either a rotation around the z-axis with frequency `ωs[i]` and duration `t`, or a rotation around the x-axis with angle `t`.
"""
function states_by_sequence(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1))
    states = SingleQubitState[s_start for _ in ωs]
    n = length(ωs)

    for pulse in pulses
        dir, t = pulse
        for i in 1:n
            rot = (dir == :z ? Rz(ωs[i]*t) : Rx(t))
            states[i] = rot * states[i]
        end
    end

    return states
end

"""
    function overlap_by_sequence(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1), s_target::SingleQubitState=SingleQubitState(1,0))

Calculates the overlap of the states obtained from `states_by_sequence` with a target state `s_target`. The overlap is computed as the squared absolute value of the inner product between each state and the target state.
"""
function overlap_by_sequence(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1), s_target::SingleQubitState=SingleQubitState(1,0))
    states = states_by_sequence(ωs, pulses, s_start=s_start) 
    overlap = [abs2(dot(state, s_target)) for state in states] 
    return overlap
end

"""
    function states_by_sequence_threaded(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1))

Calculates the states of a single qubit after applying a sequence of pulses defined by `pulses` to the initial state `s_start`, using multithreading for performance. Each pulse is either a rotation around the z-axis with frequency `ωs[i]` and duration `t`, or a rotation around the x-axis with angle `t`.
"""
function states_by_sequence_threaded(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1))
    states = SingleQubitState[s_start for _ in ωs]
    n = length(ωs)
    
    @threads for i in 1:n
        for pulse in pulses
            dir, t = pulse
            rot = (dir == :z ? Rz(ωs[i]*t) : Rx(t))
            states[i] = rot * states[i]
        end
    end

    return states
end

"""
    function overlap_by_sequence_threaded(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1), s_target::SingleQubitState=SingleQubitState(1,0))

Calculates the overlap of the states obtained from `states_by_sequence_threaded` with a target state `s_target`, using multithreading for performance. The overlap is computed as the squared absolute value of the inner product between each state and the target state.
"""
function overlap_by_sequence_threaded(ωs::Vector{BigFloat}, pulses::Vector{Tuple{Symbol, BigFloat}}; s_start::SingleQubitState=SingleQubitState(1,1), s_target::SingleQubitState=SingleQubitState(1,0))
    states = states_by_sequence_threaded(ωs, pulses, s_start=s_start)
    overlap = zeros(BigFloat, length(ωs)) 
    ThreadsX.map!(s -> abs2(dot(s, s_target)), overlap, states)
    return overlap
end

# HELPERS
"""
    function vals2sequence(vals::Vector{BigFloat})::Vector{Tuple{Symbol, BigFloat}}

Converts a vector of `BigFloat` values into a sequence of tuples, alternating between :z and :x directions. The length of the input vector must be even.
"""
function vals2sequence(vals::Vector{BigFloat})::Vector{Tuple{Symbol, BigFloat}}
    if isodd(length(vals)) 
        error("Odd sequence length!")
    end
    return [((isodd(i) ? :z : :x) , v) for (i,v) in enumerate(vals)]
end

"""
    function seq2vals(sequence::Vector{Tuple{Symbol, BigFloat}})::Vector{BigFloat}

Converts a sequence of tuples (direction, value) into a vector of `BigFloat` values, extracting only the values from the tuples.
"""
function seq2vals(sequence::Vector{Tuple{Symbol, BigFloat}})::Vector{BigFloat}
    return [t for (dir,t) in sequence]
end

end
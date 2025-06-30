# setup plotting
function setup_blochplot(;size=(500,500), fig=false, eyepos = [3,3,2], _unusedkwargs...)
    if fig == false
        fig = Figure(size = size)
    end
    ax = LScene(fig[1, 1], show_axis=true, scenekw=(limits=Rect(-1,-1,-1,2,2,2),)) 
    cam = cameracontrols(ax.scene)
    #cam.lookat[] = [0,0,0]
    cam.eyeposition[] = eyepos
    update_cam!(ax.scene, cam)
    #=
    ax = Axis3(fig[1, 1])
    xlims!(ax, (-1,1))
    ylims!(ax, (-1,1)) 
    zlims!(ax, (-1,1)) 
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.zlabel = "z"
    ax.aspect = (1.,1.,1.)
    =#
    
    # Sphere
    mesh!(ax, Sphere(Point3f(0), 1f0), color = :lightblue, alpha=0.05)
    lines!(ax, [Point3f(cos(t), sin(t), 0) for t in LinRange(0, 2pi, 50)], linewidth=0.5, color=:black, linestyle=:dash)
    lines!(ax, [Point3f(cos(t), 0, sin(t)) for t in LinRange(0, 2pi, 50)], linewidth=0.5, color=:black, linestyle=:dash)
    lines!(ax, [Point3f(0, cos(t), sin(t)) for t in LinRange(0, 2pi, 50)], linewidth=0.5, color=:black, linestyle=:dash)
    lines!(ax, [Point3f(-1.5,0,0), Point3f(+1.5,0,0)], linewidth=0.5, color=:gray, linestyle=:solid)
    lines!(ax, [Point3f(0,-1.5,0), Point3f(0,+1.5,0)], linewidth=0.5, color=:gray, linestyle=:solid)
    lines!(ax, [Point3f(0,0,-1.5), Point3f(0,0,+1.5)], linewidth=0.5, color=:gray, linestyle=:solid)

    # state labels
    text!(ax, text=["∣0⟩"], [Point3f(0,0.,1.1)], align=(-0.25, 0))
    text!(ax, text=["∣1⟩"], [Point3f(0,0.,-1.25)], align=(-0.25, 1))
    text!(ax, text=["∣+⟩"], [Point3f(1.25,0.,0)], align=(1, 0))
    text!(ax, text=["|+i⟩"], [Point3f(0,1.25,0)], align=(0, 0))
    
    return fig, ax
end

function setup_doubleplot(;size=(1500, 1200))
    fig = Figure(size=size)
    _, ax_bloch = setup_blochplot(fig=fig[1,1:2])
    ax_stereo = PolarAxis(fig[1,3])
    ax_stereo.title = "Stereographic (Northpole)"
    rlims!(ax_stereo, 0,1.05)

    return fig, ax_bloch, ax_stereo
end

function setup_tripleplot(;size=(1500, 1200))
    fig = Figure(size=size)
    _, ax_bloch = setup_blochplot(fig=fig[1:2,1:2])
    ax_stereo = PolarAxis(fig[1,3])
    ax_stereo_s = PolarAxis(fig[2,3])
    ax_stereo.title = "Stereographic (Northpole)"
    ax_stereo_s.title = "Stereographic (Southpole)"
    rlims!(ax_stereo, 0,1.05)
    rlims!(ax_stereo_s, 0,1.05)

    return fig, ax_bloch, ax_stereo, ax_stereo_s
end

function plot_axis_arrow(ax, P_rot; linewidth=0.01, color=:blue)
    arrows!(ax, [0], [0], [0],
        [P_rot[1]], [P_rot[2]], [P_rot[3]], 
        linewidth = linewidth, 
        arrowsize = Vec3f(0.05, 0.05, 0.04),
        color=color)
end

# path interpolation for a single qubit state
"""
    interpolate_path(state::SingleQubitState, ϕ::Real, Rot::Function, N::Int=50)

Interpolates a path of `N` states starting from `state` by applying the single qubit operator `Rot` with parameter `ϕ/N` for `N` steps.
The function returns a vector of `SingleQubitState` objects representing the interpolated states.
"""
function interpolate_path(state::SingleQubitState, ϕ::Real, Rot::Function, N::Int=50)
    dt = 1/N
    states = SingleQubitState[state]
    rot = Rot(ϕ * dt)
    for i in 1:N
        push!(states, rot * states[end])
    end
    return states[2:end]
end

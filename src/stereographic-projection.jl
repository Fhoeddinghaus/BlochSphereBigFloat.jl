# MOVE to BlochSphere.jl
mutable struct StereographicCoordinates{ST <: Union{PolarCoordinates, PlanarCoordinates, Complex}} <: Coordinates{BigFloat, 2}
    data::ST   # stores the actual coordinates as Polar, Planar Coordinates or as a Complex number
    sgnZ::Int8 # stores the sign (hemisphere) of the original point, independent of the pole of projection
end

function Base.length(p::StereographicCoordinates{ST}) where {ST <: Union{PolarCoordinates, PlanarCoordinates, Complex}}
    if ST == Complex
        return 2
    else
        return 3
    end
end

function Base.iterate(p::StereographicCoordinates{ST}, state::Int64=1) where {ST <: Union{PolarCoordinates, PlanarCoordinates, Complex}}
    if state > length(p)
        return nothing
    end
    if ST == Complex
        return (( state == 1 ? p.data : p.sgnZ), state+1)
    else
        if state == 3
            return (p.sgnZ, state+1)
        else
            return iterate(p.data, state)
        end
    end
end

function Base.propertynames(p::StereographicCoordinates)
    return vcat([:data, :sgnZ], [propertynames(p.data)...]) 
end

function Base.getproperty(p::StereographicCoordinates, i::Symbol)
    if i in (:data, :sgnZ)
        return getfield(p, i)
    elseif hasproperty(p.data, i)
        return getproperty(p.data, i)
    else
        error("type $(typeof(p)) has no field $i")
    end
end

function Base.convert(::Type{StereographicCoordinates{PolarCoordinates}}, s::StereographicCoordinates{ST}) where ST
    if ST == PolarCoordinates
        return s
    elseif ST == Complex # := PlanarCoordinates
        x,y = s.re, s.im
        pl = PlanarCoordinates([x,y])
    else # PlanarCoordinates
        pl = s.data
    end
    pol = convert(PolarCoordinates, pl)
    return StereographicCoordinates{PolarCoordinates}(pol, s.sgnZ)
end

function Base.convert(::Type{StereographicCoordinates{PlanarCoordinates}}, s::StereographicCoordinates{ST}) where ST
    if ST == PlanarCoordinates
        return s
    elseif ST == Complex # := PlanarCoordinates
        x,y = s.re, s.im
        pl = PlanarCoordinates([x,y])
    else # PolarCoordinates
        pl = convert(PlanarCoordinates{BigFloat}, s.data)
    end
    return StereographicCoordinates{PlanarCoordinates}(pl, s.sgnZ)
end

function Base.convert(::Type{StereographicCoordinates{Complex}}, s::StereographicCoordinates{ST}) where ST
    if ST == Complex
        return s
    elseif ST == PlanarCoordinates
        z = s.x + im * s.y
    else # PolarCoordinates
        pl = convert(PlanarCoordinates{BigFloat}, s.data)
        z = pl.x + im * pl.y
    end
    return StereographicCoordinates{Complex}(z, s.sgnZ)
end

# ↦ stereo
function projection(::Type{StereographicCoordinates}, point::SphereCoordinates; proj=:north)
    # get hemisphere
    sgnZ = sign(π/2-point.polar)
    #(0 <= point.polar < π/2 ? +1 : (point.polar == π/2 ? 0 : -1))
    
    # special case:
    if proj in (:both, :abs)
        proj = (sgnZ > 0 ? :south : :north)
    end
    
    if proj in (:top, :north, :topdown, :n, :t, +1)
        R = cot(point.polar/2)
    elseif proj in (:bottom, :south, :bottomup, :s, :b, -1)
        R = cot((π-point.polar)/2)
    else
        error("proj has to be one of (:top, :north, :topdown, :n, :t, +1) or (:bottom, :south, :bottomup, :s, :b, -1) or (:both, :abs)")
    end
    p = PolarCoordinates(R, point.azi)
    return StereographicCoordinates{PolarCoordinates}(p, sgnZ)
end

function projection(::Type{StereographicCoordinates}, point::SpatialCoordinates; proj=:north)
    # get hemisphere
    sgnZ = sign(point.z)
    
    if proj in (:top, :north, :topdown, :n, :t, +1)
        sz = point.z
    elseif proj in (:bottom, :south, :bottomup, :s, :b, -1)
        sz = -point.z
    elseif proj in (:both, :abs)
        sz = -abs(point.z)
    else
        error("proj has to be one of (:top, :north, :topdown, :n, :t, +1) or (:bottom, :south, :bottomup, :s, :b, -1) or (:both, :abs)")
    end
    p = PlanarCoordinates([
        point.x / (1 - sz),
        point.y / (1- sz)
    ])
    return StereographicCoordinates{PlanarCoordinates}(p, sgnZ)
end

function projection(::Type{StereographicCoordinates}, point::SingleQubitState; proj=:north)
    # get hemisphere
    sgnZ = sign(round(abs2(point[1]) - abs2(point[2]), digits=numdigits()-1)) 
    
    # special case:
    if proj in (:both, :abs)
        proj = (sgnZ > 0 ? :south : :north)
    end
    
    if proj in (:top, :north, :topdown, :n, :t, +1)
        # ref: https://arxiv.org/abs/quant-ph/0201014v2, Wiki
        # c = exp(im*φ) * cot(θ/2) = a/b; north is defined as conj(c=a/b)
        c = conj(point[1]/point[2])
    elseif proj in (:bottom, :south, :bottomup, :s, :b, -1)
        # equivalent to switching |0⟩, |1⟩. south is defined as c=b/a
        c = point[2]/point[1]
    else
        error("proj has to be one of (:top, :north, :topdown, :n, :t, +1) or (:bottom, :south, :bottomup, :s, :b, -1) or (:both, :abs)")
    end
    if isnan(c)
        c = Inf + Inf*im
    end
    return StereographicCoordinates{Complex}(c, sgnZ)
end

# stereo ↦
function projection(::Type{SphereCoordinates}, s::StereographicCoordinates)
    s = convert(StereographicCoordinates{PolarCoordinates}, s)
    R, azi = s.r, s.azi
    polar = 2 * atan(1/R)
    sgnZ = (R <= 1 ? +1 : -1) * s.sgnZ
    polar = (sgnZ < 0 ? polar : π - polar)
    return SphereCoordinates(polar, azi)
end

function projection(::Type{SpatialCoordinates}, s::StereographicCoordinates)
    s = convert(StereographicCoordinates{PlanarCoordinates}, s)
    X = s.data
    sgnZ = (norm(X) <= 1 ? -1 : +1) * s.sgnZ
    
    rr = X[1]^2 + X[2]^2 + 1 # rr - 2 = X² + Y² - 1
    x = [2 * X[1] / rr,
        2 * X[2] / rr,
        s * (rr-2) / rr]
    return SpatialCoordinates(x)
end

function projection(::Type{SingleQubitState}, s::StereographicCoordinates)
    # c = a/b
    # c = exp(im*φ) * cot(θ/2) = phase * radius = a/b
    # c/exp(im*φ) = cot(θ/2)  = R == sqrt(c*c')
    # θ = 2 * atan(1/R) == 2 * atan(1/sqrt(c*c'))
    # φ = imag(log(c/sqrt(c*c')))
    # phase = c/sqrt(c*c')
    
    # ⟹ !!!
    # a = cos(θ/2) == cos(atan(1/sqrt(c*c'))) 
    #   = 1/sqrt(1 + 1/(c*c'))
    
    # b = phase * sin(θ/2)
    #   = 1/(sqrt(1 + 1/(c*c')) * c')
    # for bottom: use conj(c) and switch a,b
    s = convert(StereographicCoordinates{Complex}, s)
    c = s.data
    
    # |c|^sgnZ: for :bottom ∈ [0,1), :both = 1, :top ∈ (1,∞)
    directed_R = norm(c)^(s.sgnZ)
    if directed_R < 1
        # bottom
        c = conj(c)
        a = 1/(sqrt(1 + 1/(c*c')) * c')
        b = 1/sqrt(1 + 1/(c*c')) 
    else
        # top
        a = 1/sqrt(1 + 1/(c*c')) 
        b = 1/(sqrt(1 + 1/(c*c')) * c')
    end
    return SingleQubitState([a,b])
end
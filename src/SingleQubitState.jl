mutable struct SingleQubitState <: StaticArrays.FieldVector{2, Complex{BigFloat}}
    # defines |ψ⟩ = a|0⟩ + b|1⟩ with standard computational basis. ψ = (a,b)ᵀ
    a::Complex{BigFloat}
    b::Complex{BigFloat}
    
    function SingleQubitState(a::Complex{BigFloat}, b::Complex{BigFloat})
        # 1. normalize
        N = sqrt(a*a' + b*b')
        if N != 1
            a,b = [a,b]/N
        end
        
        # 2. extract global phase s.t. Im(a) == 0
        # first, try setting b.im to 0, then a.im to 0 for more conventional notations
        
        if b.im != 0
            a,b = round.(exp(-im * log(b).im) * [a,b], digits=numdigits()-1)
        end
        # for form |ψ⟩ = cos(θ/2)|0⟩ + exp(im φ) sin(θ/2)|1⟩
        if a.im != 0
            a,b = round.(exp(-im * log(a).im) * [a,b], digits=numdigits()-1)
        end
        new(a,b)
    end
    SingleQubitState(a::ComplexF64, b::ComplexF64) = SingleQubitState(Complex{BigFloat}(a),Complex{BigFloat}(b))
    SingleQubitState(a::T, b::T) where {T <: Real} = SingleQubitState(Complex{BigFloat}(a),Complex{BigFloat}(b))
    SingleQubitState(s::Vector{T}) where {T <: Real} = SingleQubitState(s[1], s[2])
    SingleQubitState(s::Vector{Complex{T}}) where {T <: Real} = SingleQubitState(s[1], s[2])
end

Base.show(io::IO, s::SingleQubitState) = print_type_styled(io, 
    "SingleQubitState", 
    "{BigFloat}", 
    "(a=$(s.a), b=$(s.b))"
)

function BlochVec(s::SingleQubitState)::SpatialCoordinates{BigFloat}
    x = 2 * real(s[1] * s[2]')
    y = -2 * imag(s[1] * s[2]')
    #z = real(s[1] * s[1]' - s[2] * s[2]')
    z = abs2(s[1]) - abs2(s[2])
    return SpatialCoordinates([x,y,z])
end
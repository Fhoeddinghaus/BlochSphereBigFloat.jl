p1 = [0 1; 
    1 0]

PauliX = p1

p2 = [0  -im; 
    im  0]
PauliY = p2

p3 = [1  0; 
    0 -1]
PauliZ = p3

P = [p1, p2, p3]

"""
    ep(x::Vector) = exp(-im/2 * sum(x * P)) 

Calculates the matrix exponential `e^(-im/2 * ∑ xᵢ⋅Pauliᵢ)` for a vector `x` where `Pauli[X,Y,Z]` are the Pauli matrices.
"""
ep(x) = exp(-im/2 * sum(x * P))

"""
    Rotation operator on ℂ² around x-axis (in Bloch picture) by angle `β`.
"""
Rx(β) = cos(β/2) * I - im * sin(β/2) * PauliX
# more precise than
#Rx(β) = ep(β * [1 0 0])

"""
    Rotation operator on ℂ² around y-axis (in Bloch picture) by angle `β`.
"""
Ry(β) = cos(β/2) * I - im * sin(β/2) * PauliY
# more precise than
#Ry(β) = ep(β * [0 1 0])

"""
    Rotation operator on ℂ² around z-axis (in Bloch picture) by angle `γ`.
"""
Rz(γ) = cos(γ/2) * I - im * sin(γ/2) * PauliZ
# more precise than
#Rz(γ) = ep(γ * [0 0 1])

"""
    numdigits()

Returns the number of digits of the current precision for `BigFloat` numbers, which is determined by the machine epsilon.
The precision can be adjusted using `setprecision(BigFloat, num_bits)` for `num_bits` bits.
"""
numdigits() = max(0, ceil(Int, -log10(eps(BigFloat))))
module SingleSpin

"""
    X, Y, Z

Single-spin Pauli matrices.
"""
X = [0 1; 1 0] ./ 2
Y = [0 -im; im 0] ./ 2
Z = [1 0; 0 -1] ./ 2

"""
    h_pulse(ω1, φ)

Generate a Hamiltonian for a pulse with a given angular frequency (ω1 = B1/γ)
and pulse phase.

The angular frequency should be specified in rad/s. The pulse phase can be
either specified as an angle in radians, or as one of the following symbols:
`:x`, `:y`, `:_x`, or `:_y`, which correspond respectively to 0, π/2, π, and
3π/2.
"""
function h_pulse(ω1, φ)
    if φ === :x
        φ = 0
    elseif φ === :y
        φ = pi/2
    elseif φ === :_x
        φ = pi
    elseif φ === :_y
        φ = 3pi/2
    end

    return ω1 * (cos(φ) * X + sin(φ) * Y)
end

"""
    h_offset(Ω)

Generate a Hamiltonian for a spin offset.

The offset frequency should be specified in Hz.
"""
function h_offset(Ω)
    Ω * 2 * pi * Z
end

end  # module SingleSpin

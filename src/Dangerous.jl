module Dangerous

import LinearAlgebra as LA
import FFTW

function propagate(rho, H, t)
    U = exp(-im * H * t)
    U * rho * U'
end


module SingleSpin

"""
    X, Y, Z

Single-spin Pauli matrices.
"""
X = [0 1; 1 0] ./ 2
Y = [0 -im; im 0] ./ 2
Z = [1 0; 0 -1] ./ 2

"""
    h_pulse(ω1, θ)

Generate a Hamiltonian for a pulse with a given angular frequency (ω1 = B1/γ)
and pulse axis.

The angular frequency should be specified in rad/s. The pulse axis can be
either specified as an angle in radians, or as one of the following symbols:
`:x`, `:y`, `:_x`, or `:_y`, which correspond respectively to 0, π/2, π, and
3π/2.
"""
function h_pulse(ω1, θ)
    if θ === :x
        θ = 0
    elseif θ === :y
        θ = pi/2
    elseif θ === :_x
        θ = pi
    elseif θ === :_y
        θ = 3pi/2
    end

    return ω1 * (cos(θ) * X + sin(θ) * Y)
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

function detect(rho, h_free, dt, npoints)
    fid::Array{ComplexF64} = zeros(npoints)
    for i in 1:npoints
        fid[i] = -LA.tr(rho * SingleSpin.X) - im * LA.tr(rho * SingleSpin.Y)
        old_norm = LA.norm(rho)
        rho = propagate(rho, h_free, dt)
        new_norm = LA.norm(rho)
        if abs(old_norm - new_norm) >= 1e-14
            error("Norm not preserved during propagation step: $old_norm -> $new_norm")
        end
    end

    # Fake relaxation
    k = 0.01
    window = exp.(-k .* (1:npoints))
    fid = fid .* window

    fid
end

# Get the spectrum of a single spin with chemical shift cs
# TODO: How do you document Julia functions?
function spectrum_single_spin(sfo1, cs, aq, td)
    # Start with a 90(-y) pulse
    rho = SingleSpin.Z
    h_free = SingleSpin.h_offset(sfo1 * cs)
    rho = propagate(rho, SingleSpin.h_pulse(pi/2, :_y), 1)

    # Then detect and Fourier transform
    dw = aq / td
    fid = detect(rho, h_free, dw, td)
    x = FFTW.fftshift(FFTW.fftfreq(td, 1/dw)) ./ sfo1
    y = FFTW.fftshift(FFTW.fft(fid))
    x, y
end

end # module Dangerous

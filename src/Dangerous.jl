module Dangerous

import LinearAlgebra as LA
import FFTW

# Single spin Pauli matrices

X = [0 1; 1 0] ./ 2
Y = [0 -im; im 0] ./ 2
Z = [1 0; 0 -1] ./ 2

function propagate(rho, H, t)
    U = exp(-im * H * t)
    U * rho * U'
end

function h_pulse(flip_angle, axis)
    axis = lowercase(axis)
    if axis == "x"
        flip_angle * X
    elseif axis == "y"
        flip_angle * Y
    elseif axis == "z"
        flip_angle * Z
    elseif axis == "-x"
        -flip_angle * X
    elseif axis == "-y"
        -flip_angle * Y
    elseif axis == "-z"
        -flip_angle * Z
    else
        error("unknown axis '$axis'")
    end
end

function h_offset(omega)
    omega * 2 * pi * Z
end

function detect(rho, h_free, dt, npoints)
    fid::Array{ComplexF64} = zeros(npoints)
    for i in 1:npoints
        fid[i] = -LA.tr(rho * X) - im * LA.tr(rho * Y)
        old_norm = LA.norm(rho)
        rho = propagate(rho, h_free, dt)
        new_norm = LA.norm(rho)
        if abs(old_norm - new_norm) >= 1e-3
            error("Norm not preserved: $old_norm -> $new_norm")
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
    rho = Z
    h_free = h_offset(sfo1 * cs)
    rho = propagate(rho, h_pulse(pi/2, "-y"), 1)

    # Then detect and Fourier transform
    dw = aq / td
    fid = detect(rho, h_free, dw, td)
    x = FFTW.fftshift(FFTW.fftfreq(td, 1/dw)) ./ sfo1
    y = FFTW.fftshift(FFTW.fft(fid))
    x, y
end

export X, Y, Z, propagate

end # module Dangerous

module Dangerous

import LinearAlgebra as LA
import FFTW

include("SingleSpin.jl")


function propagate(rho, H, t)
    U = exp(-im * H * t)
    U * rho * U'
end

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

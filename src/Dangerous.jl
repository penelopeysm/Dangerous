module Dangerous

import LinearAlgebra as LA
import FFTW

include("Units.jl")
include("Nuclei.jl")
include("SingleSpin.jl")

module SpinSystem

# struct System
#     nuclei::Array{SingleSpin.Nucleus}
# end

end # module SpinSystem


function propagate(ρ, H, t)
    # TODO: Is there a way to accumulate the propagators so that we can halve
    # the number of matrix multiplications?
    U = exp(-im * H * t)
    return U * ρ * U'
end

function detect(ρ, Hfree, dwell, npoints)
    fid::Array{ComplexF64} = zeros(npoints)
    for i in 1:npoints
        fid[i] = -LA.tr(ρ * SingleSpin.X) - im * LA.tr(ρ * SingleSpin.Y)
        old_norm = LA.norm(ρ)
        ρ = propagate(ρ, Hfree, dwell)
        new_norm = LA.norm(ρ)
        if abs(old_norm - new_norm) >= 1e-14
            error("Norm not preserved during propagation step: $old_norm -> $new_norm")
        end
    end

    # Fake relaxation
    k = 0.01
    window = exp.(-k .* (1:npoints))
    fid = fid .* window

    return fid
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
    return x, y
end

end # module Dangerous

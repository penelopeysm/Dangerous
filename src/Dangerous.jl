module Dangerous

import LinearAlgebra as LA
import FFTW

include("units.jl")
include("nuclei.jl")
include("spin_system.jl")
include("hamiltonian.jl")


function propagate(ρ, H, t)
    # TODO: Is there a way to accumulate the propagators so that we can halve
    # the number of matrix multiplications?
    U = exp(-im * H * t)
    return U * ρ * U'
end

function detect(sys, nuc, ρ, dwell, npoints)
    fid = zeros(ComplexF64, npoints)
    x, y = Hamiltonian.detection_operators(sys, nuc)
    Hfree = Hamiltonian.h_free(sys)
    for i in 1:npoints
        fid[i] = -LA.tr(ρ * x) - im * LA.tr(ρ * y)
        old_norm = LA.norm(ρ)
        ρ = propagate(ρ, Hfree, dwell)
        new_norm = LA.norm(ρ)
        if abs(old_norm - new_norm) >= 1e-14
            error("Norm not preserved during propagation step: $old_norm -> $new_norm")
        end
    end

    # Fake relaxation
    k = 0.005
    window = exp.(-k .* (1:npoints))
    fid = fid .* window

    return fid
end

# Get the spectrum of a spin system
function spectrum(sys, aq, td)
    # Start with a 90(-y) pulse
    ρ = Hamiltonian.ρ_eq(sys)
    ρ = propagate(ρ, Hamiltonian.h_pulse(sys, Nuclei.H1, pi/2, :_y), 1)

    # Then detect and Fourier transform
    dw = aq / td
    fid = detect(sys, Nuclei.H1, ρ, dw, td)
    x = FFTW.fftshift(FFTW.fftfreq(td, 1/dw)) / sys.magnetic_field
    y = FFTW.fftshift(FFTW.fft(fid))
    return x, y
end

end # module Dangerous

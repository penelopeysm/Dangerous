module Dangerous

import LinearAlgebra as LA
import FFTW

include("units.jl")
include("nuclei.jl")
include("spin_system.jl")
include("hamiltonian.jl")

using .Nuclei: γ
using Unitful: Time, Frequency, @u_str, ustrip, uconvert, NoUnits

function propagate(ρ, H::Matrix{<:Frequency}, t::Time)
    # TODO: Is there a way to accumulate the propagators so that we can halve
    # the number of matrix multiplications?
    U = exp(-im * H * t)
    return U * ρ * U'
end

function detect_fid(sys, nuc, ρ, dwell::Time, npoints)
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

function detect_spectrum(sys, nuc, ρ, dwell::Time, npoints)
    fid = detect_fid(sys, nuc, ρ, dwell, npoints)
    # Units are a bit of a faff here... but it's worth getting them right
    # The FFT library doesn't convert seconds into Hz correctly, so we do it
    # manually
    x_hz = FFTW.fftshift(FFTW.fftfreq(npoints, 1 / ustrip(u"s", dwell))) * u"Hz"
    # Then convert to ppm. We need to strip the units so that the plotting
    # library doesn't get confused
    x_ppm = uconvert.(NoUnits, 1e6 * x_hz / (sys.magnetic_field * γ(nuc)))
    # This one's easy!
    y = FFTW.fftshift(FFTW.fft(fid))
    return x_ppm, y
end

# Get the spectrum of a spin system
function spectrum(sys, aq::Time, td)
    # Start with a 90(-y) pulse. Note the hacky dimensions for an instantaneous pulse,
    # this should be fixed at some point.
    ρ = Hamiltonian.ρ_eq(sys)
    H_pulse = Hamiltonian.h_pulse(sys, Nuclei.H1, (π / 2)u"Hz", :_y)
    ρ = propagate(ρ, H_pulse, 1u"s")

    # Then detect and Fourier transform
    dw = aq / td
    return detect_spectrum(sys, Nuclei.H1, ρ, dw, td)
end

end # module Dangerous

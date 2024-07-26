function propagate(ρ, H::Matrix{<:Frequency}, t::Time)
    # TODO: Is there a way to accumulate the propagators so that we can halve
    # the number of matrix multiplications?
    U = exp(-im * H * t)
    return U * ρ * U'
end

function detect_fid(sys, nuc, ρ, dwell::Time, npoints)
    fid = zeros(ComplexF64, npoints)
    x, y = detection_operators(sys, nuc)
    Hfree = h_free(sys)
    for i in 1:npoints
        fid[i] = -tr(ρ * x) - im * tr(ρ * y)
        old_norm = norm(ρ)
        ρ = propagate(ρ, Hfree, dwell)
        new_norm = norm(ρ)
        if abs(old_norm - new_norm) >= 1e-10
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
    resonance_frequency = sys.magnetic_field * γ(nuc)
    println("resonance_frequency: $resonance_frequency")
    # The FFT library doesn't convert seconds into Hz correctly, so we do it
    # manually. We then need to tack on the transmitter offset
    offset_hz = sys.transmitter_offset[nuc] * resonance_frequency / 1e6
    println("offset_hz: $(uconvert(u"Hz", offset_hz))")
    x_hz = fftshift(fftfreq(npoints, 1 / ustrip(u"s", dwell))) * u"Hz" .+ offset_hz
    println("x_hz: $(minimum(collect(x_hz))) -> $(maximum(collect(x_hz)))")
    # println(collect(x_hz))
    # println(collect(x_hz / resonance_frequency))
    # Then convert to ppm. We need to strip the units so that the plotting
    # library doesn't get confused
    # Need to collect(x_hz) first otherwise it leads to a weird bug with units
    x_ppm = uconvert.(NoUnits, collect(x_hz) / resonance_frequency) * 1e6
    println("x_ppm: $(minimum(collect(x_ppm))) -> $(maximum(collect(x_ppm)))")
    # This one's easy!
    y = fftshift(fft(fid))
    return x_ppm, y
end

"""
    zg(sys, nuc, sw, td)

1D pulse–acquire experiment.

# Arguments

- `sys::`[`Dangerous.System`](@ref): the spin system to simulate 
- `nuc::`[`Dangerous.Nucleus`](@ref): the nucleus to detect
- `sw`: the spectral window (in ppm)
- `td`: the number of points in the FID
"""
function zg(sys, nuc, sw, td)
    # Start with a 90(-y) pulse. Note the hacky dimensions for an instantaneous pulse,
    # this should be fixed at some point.
    ρ = ρ_eq(sys)
    H_pulse = h_pulse(sys, nuc, (π / 2)u"Hz", :_y)
    ρ = propagate(ρ, H_pulse, 1u"s")

    # Then detect and Fourier transform
    dw = 1 / ((sw * 1e-6) * sys.magnetic_field * γ(nuc))
    println("aq: $(uconvert(u"s", dw * td))")
    return detect_spectrum(sys, nuc, ρ, dw, td)
end

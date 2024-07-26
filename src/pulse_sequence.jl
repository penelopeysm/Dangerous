"""
    SpinSystemState(ρ_0, U)

This struct holds the state of a spin system, but without explicitly
propagating the density matrix. Instead, the initial density matrix and the
accumulated propagator are separately stored. Only at the end of the pulse
sequence is the propagation performed. This allows us to halve the number of
matrix multiplications being done.
"""
struct SpinSystemState
    ρ_0
    U
end

"""
    SpinSystemState(ρ_0)

Create a `SpinSystemState` with an initial density matrix `ρ_0` and no
propagation yet (i.e. `U` is the identity matrix).
"""
function SpinSystemState(ρ_0)
    n = size(ρ_0, 1)
    SpinSystemState(ρ_0, Matrix{ComplexF64}(I, n, n))
end

"""
    propagate(ρ, H, t)

Propagate the density matrix `ρ` under the Hamiltonian `H` for a time `t`.
"""
function propagate(ρ, H::Matrix{<:Frequency}, t::Time)
    U = exp(-im * H * t)
    return U * ρ * U'
end

"""
    propagate(state::SpinSystemState, H, t)

Propagate the uncollapsed spin system state under the Hamiltonian `H` for a
time `t`.
"""
function propagate(state::SpinSystemState, H::Matrix{<:Frequency}, t::Time)
    this_U = exp(-im * H * t)
    return SpinSystemState(state.ρ_0, this_U * state.U)
end

"""
    collapse(state)

Collapse the intermediate spin system state by actually performing the
propagation on the density matrix.
"""
function collapse(state::SpinSystemState)
    return state.U * state.ρ_0 * state.U'
end

"""
    detect_fid(sys, nuc, ρ, dwell, npoints)

Detect a free induction decay.
"""
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
    T2 = 0.5u"s"
    window = @. exp((1:npoints) * -dwell / T2)
    fid = fid .* window

    return fid
end

"""
    detect_spectrum(sys, nuc, ρ, dwell, npoints)

Detect a spectrum. Returns a tuple of the frequency axis and the spectrum.
"""
function detect_spectrum(sys, nuc, ρ, dwell::Time, npoints)
    fid = detect_fid(sys, nuc, ρ, dwell, npoints)
    # Units are a bit of a faff here... but it's worth getting them right
    resonance_frequency = sys.magnetic_field * γ(nuc)
    # The FFT library doesn't convert seconds into Hz correctly, so we do it
    # manually. We then need to tack on the transmitter offset
    offset_hz = sys.transmitter_offset[nuc] * resonance_frequency / 1e6
    x_hz = fftshift(fftfreq(npoints, 1 / ustrip(u"s", dwell))) * u"Hz" .+ offset_hz
    # println(collect(x_hz))
    # println(collect(x_hz / resonance_frequency))
    # Then convert to ppm. We need to strip the units so that the plotting
    # library doesn't get confused
    # Need to collect(x_hz) first otherwise it leads to a weird bug with units
    x_ppm = uconvert.(NoUnits, collect(x_hz) / resonance_frequency) * 1e6
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
    ρ = SpinSystemState(ρ_eq(sys))
    # Start with a 90(-y) pulse. Note the hacky dimensions for an instantaneous pulse,
    # this should be fixed at some point.
    ρ = propagate(ρ, h_pulse(sys, nuc, (π / 2)u"Hz", :_y), 1u"s")
    ρ = collapse(ρ)

    # Then detect and Fourier transform
    dw = 1 / ((sw * 1e-6) * sys.magnetic_field * γ(nuc))
    @info "zg: aq = $(uconvert(u"s", dw * td))"
    @info "zg: dw = $(uconvert(u"µs", dw))"
    return detect_spectrum(sys, nuc, ρ, dw, td)
end

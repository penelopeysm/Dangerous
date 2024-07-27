using Dangerous
using Logging
using Unitful: @u_str

sys = System(
    14.1u"T",
    Dict(H1 => 4.7),
    [H1, H1],
    [1.5, 9],
    [0 30; 0 0]u"Hz"
)

# Either use the builtin functions; this is a simple spin echo

@pulse_sequence sys begin
    pulse_instant(H1, π / 2, :_y)
    delay(5u"ms")
    pulse_instant(H1, π, :x)
    delay(5u"ms")
    spec = detect_spectrum_1d(H1, 15, 32768)
end

using GLMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity", title="Spectrum")
ax.xreversed = true
lines!(ax, spec.x, real.(spec.y), color=:blue, label="Real")
# lines!(ax, spec.x, imag.(spec.y), color=:orange, label="Imaginary")
axislegend()
display(f)

# You can define your own functions if you want, but there'll be LSP errors

if false  # Disabled, so that this file can be imported easily
    function my_pulse_instant(ρ, sys, nuc, angle, phase)
        return Dangerous.propagate(ρ, Dangerous.h_pulse(sys, nuc, angle * u"Hz", phase), 1u"s")
    end

    @pulse_sequence sys begin
        my_pulse_instant(H1, π / 2, :_y)
        spec = detect_spectrum_1d(H1, 15, 32768)
    end
end

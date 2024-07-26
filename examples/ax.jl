using Dangerous: Nuclei, SpinSystem, spectrum
using Unitful: @u_str

sys = SpinSystem.System(
    14.1u"T",
    Dict(Nuclei.H1 => 4.7),
    [Nuclei.H1, Nuclei.H1],
    [1.5, 7],
    [0 30; 0 0]u"Hz"
)

freq, spec = spectrum(sys, 20, 32768)

using CairoMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity", title="Spectrum")
ax.xreversed = true
lines!(ax, freq, real.(spec), color=:blue, label="Real")
lines!(ax, freq, imag.(spec), color=:orange, label="Imaginary")
axislegend()
display(f)

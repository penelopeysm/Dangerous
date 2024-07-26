using Dangerous: Nuclei, SpinSystem, spectrum
using Unitful: @u_str

sys = SpinSystem.System(
    22.1u"T",
    [Nuclei.H1],
    [1.5],
    [0.0u"Hz";;]
)

freq, spec = spectrum(sys, 3u"s", 32768)

using CairoMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity", title="Spectrum")
ax.xreversed = true
lines!(ax, freq, real.(spec), color=:blue, label="Real")
lines!(ax, freq, imag.(spec), color=:orange, label="Imaginary")
axislegend()
display(f)

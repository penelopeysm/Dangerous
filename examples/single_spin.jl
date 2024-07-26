using Dangerous: Nuclei, SpinSystem, spectrum

sys = SpinSystem.System(
    600,
    [Nuclei.H1],
    [1000],
    [0;;]
)

freq, spec = spectrum(sys, 3, 32768)

using CairoMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity", title="Spectrum")
ax.xreversed = true
lines!(ax, freq, real.(spec), color=:blue, label="Real")
lines!(ax, freq, imag.(spec), color=:orange, label="Imaginary")
axislegend()
display(f)

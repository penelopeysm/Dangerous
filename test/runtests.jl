import Dangerous as D

freq, spec = D.spectrum_single_spin(
    600, -1.5, 3, 32768
)

using GLMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity",
          title="Spectrum of a single spin")
ax.xreversed = true
lines!(ax, freq, real.(spec), color=:blue, label="Real")
lines!(ax, freq, imag.(spec), color=:orange, label="Imaginary")
axislegend()
display(f)

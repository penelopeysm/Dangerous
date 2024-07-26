using Dangerous
using Unitful: @u_str

sys = System(
    14.1u"T",   # 600 MHz 1H
    Dict(H1 => 1.7),
    [H1],
    [1.5],
    [0.0u"Hz";;]
)

freq, spec = zg(sys, H1, 20, 32768)

using CairoMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity", title="Spectrum")
ax.xreversed = true
lines!(ax, freq, real.(spec), color=:blue, label="Real")
lines!(ax, freq, imag.(spec), color=:orange, label="Imaginary")
axislegend()
display(f)

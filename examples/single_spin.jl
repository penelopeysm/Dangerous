using Dangerous
using Logging
using Unitful: @u_str

global_logger(ConsoleLogger(stderr, Info))

sys = System(
    14.1u"T",
    Dict(H1 => 1.7),
    [H1],
    [1.5],
    [0.0u"Hz";;]
)

spec = zg(sys, H1, 20, 32768)

using CairoMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity", title="Spectrum")
ax.xreversed = true
lines!(ax, spec.x, real.(spec.y), color=:blue, label="Real")
lines!(ax, spec.x, imag.(spec.y), color=:orange, label="Imaginary")
axislegend()
display(f)

# dangerous

Small Julia library to simulate 1D and 2D nuclear magnetic resonance pulse programmes.


> [!WARNING]
> This is a toy project and is specifically named to reflect that.
> Please don't take it seriously.
> In particular, please don't use it.

## Why is it dangerous?

I ran the phrase 'nuclear magnetic resonance' through Google Translate 500 times and got back 'it is dangerous'.


## What can you do with it?

If you disobey the instructions in the warning above, you can sometimes generate nice Lorentzian lineshapes.

```julia
import Dangerous as D

sys = D.SpinSystem.System(
    600,
    [D.Nuclei.H1],
    [1000],
    [0;;]
)

freq, spec = D.spectrum(
    sys, 3, 32768
)

using GLMakie
f = Figure()
ax = Axis(f[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity", title="Spectrum")
ax.xreversed = true
lines!(ax, freq, real.(spec), color=:blue, label="Real")
lines!(ax, freq, imag.(spec), color=:orange, label="Imaginary")
axislegend()
display(f)
```

![1 spin](https://github.com/user-attachments/assets/473a3535-c4fe-47d8-a2b4-6d0bd25588b2)


Alternatively, if you want a couple of spins, change the definition of the spin system:

```julia
sys = D.SpinSystem.System(
    600,
    [D.Nuclei.H1, D.Nuclei.H1],
    [1, -2],
    [0 30; 0 0]
)
```

![2 spins](https://github.com/user-attachments/assets/82a35c7a-fcf2-4c23-ad7b-d9148ad5c1f4)


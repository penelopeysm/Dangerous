# Dangerous.jl

_Dangerous_ is a Julia package for simulation of very basic nuclear magnetic resonance (NMR) experiments.

If you are trying to use this package, you probably shouldn't.
The main reason why this exists is because I needed to learn Julia.

Nevertheless, if you really want to, this documentation will hopefully teach you how to use it.

## Installation

_Dangerous_ must be installed from GitHub:

```julia
using Pkg; Pkg.add(url="https://github.com/penelopeysm/Dangerous.jl")
```

Note that Dangerous does not currently provide any functionality for plotting, it just gives you the chemical shifts and the intensities at each point (i.e. the spectrum).

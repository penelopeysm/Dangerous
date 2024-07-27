# Dangerous.jl

_Dangerous_ is a Julia package for simulation of very basic nuclear magnetic resonance (NMR) experiments.

If you are trying to use this package, you probably shouldn't.
The main reason why this exists is because I needed to learn Julia.

Nevertheless, if you really want to, this documentation will hopefully teach you how to use it.

## Installation

Dangerous must be installed from GitHub:

```julia
using Pkg; Pkg.add(url="https://github.com/penelopeysm/Dangerous.jl")
```

Note that Dangerous does not currently provide any functionality for plotting, it just gives you the chemical shifts and the intensities at each point (i.e. the spectrum).

## Basic usage

Broadly speaking, it's not too dissimilar from Spinach: you will first set up a spin system (see [`Dangerous.System`](@ref)) and then run some experiment on it (like [`Dangerous.zg`](@ref)).

I don't want to document this in too much detail because the interface is liable to change, so you can look at the `examples` folder for some examples of how to use Dangerous.

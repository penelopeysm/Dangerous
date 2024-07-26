module Dangerous

using DimensionfulAngles
using DimensionfulAngles: Periodic, radᵃ, 𝐀
using FFTW: fft, fftshift, fftfreq
using LinearAlgebra: tr, norm, I, kron
using Logging: @debug, @info, @warn, @error
using Unitful
using Unitful: 𝐈, 𝐓, 𝐌, Time, Frequency, NoUnits
using UnitfulEquivalences

export H1, C13, F19, N15, P31
export System
export zg

include("units.jl")
include("nuclei.jl")
include("spin_system.jl")
include("hamiltonian.jl")
include("pulse_sequence.jl")

end # module Dangerous

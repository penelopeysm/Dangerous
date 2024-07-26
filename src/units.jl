# NB: LSP errors on derived dimensions are spurious, see
# https://github.com/julia-vscode/StaticLint.jl/issues/381

# rad / T
@derived_dimension MagnetogyricAngular (𝐀 * 𝐈 * 𝐓 * 𝐌^-1) true
# Hz / T
@derived_dimension Magnetogyric (𝐈 * 𝐓 * 𝐌^-1) true
# T
@derived_dimension MagneticField (𝐌 * 𝐈^-1 * 𝐓^-2) true

struct PeriodicNMR <: Equivalence end
@eqrelation PeriodicNMR (MagnetogyricAngular / Magnetogyric = 2π * radᵃ)

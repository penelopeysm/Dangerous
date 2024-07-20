module Units

using Unitful
import Unitful: 𝐈, 𝐓, 𝐌
using UnitfulEquivalences
import DimensionfulAngles: Periodic, radᵃ, 𝐀

# LSP errors on these macros are spurious and can be safely ignored
@derived_dimension MagnetogyricAngular (𝐀*𝐈*𝐓*𝐌^-1) true
@derived_dimension Magnetogyric (𝐈*𝐓*𝐌^-1) true

struct PeriodicNMR <: Equivalence end
@eqrelation PeriodicNMR (MagnetogyricAngular/Magnetogyric = 2π * radᵃ)

export PeriodicNMR

end # module Units

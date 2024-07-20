module Nuclei

using ..Units

using Unitful
import UnitfulEquivalences: uconvert
using DimensionfulAngles

@enum Nucleus begin
    H1
    C13
    F19
    N15
    P31
end

function γ(n::Nucleus)
    from_bruker = if n == H1
        26.7522208e7u"radᵃ/T/s"
    elseif n == C13
        6.728286e7u"radᵃ/T/s"
    elseif n == F19
        25.16233e7u"radᵃ/T/s"
    elseif n == N15
        -2.7126189e7u"radᵃ/T/s"
    elseif n == P31
        10.8384e7u"radᵃ/T/s"
    else
        error("Unknown nucleus: '$n'")
    end
    return uconvert(u"MHz/T", from_bruker, PeriodicNMR())
end

export γ

end # module Nuclei

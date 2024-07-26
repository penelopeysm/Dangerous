module SpinSystem

using LinearAlgebra: diag
using Unitful: Frequency, ustrip
using ..Nuclei: Nucleus
using ..Units: MagneticField

"""
    System(
        magnetic_field::MagneticField,
        nuclei::Vector{Nucleus},
        chemical_shifts::Vector{Float64},
        couplings::Matrix{Unitful.Frequency}
    )

A complete specification of a spin system.

TODO: The magnetic field is not part of the physical spin system, but rather
the spectrometer / experiment setup. It should logically be moved elsewhere.

Fields:
- `magnetic_field::MagneticField`: The magnetic field strength, e.g. `22.1u"T"`.
- `nuclei::Vector{Nucleus}`: The nuclei in the system.
- `chemical_shifts::Vector{Float64}`: The chemical shifts of the nuclei in the system, in ppm. The length of this array must match the length of `nuclei`.
- `couplings::Matrix{Float64}`: The coupling constants between the nuclei in the system. The size of this matrix must be `(length(nuclei), length(nuclei))`. Note that the diagonal elements of this matrix must be zero, and for each pair of off-diagonal elements (i, j) and (j, i), one of them must be zero. The easiest way to generate this matrix is to initialise it with `zeros(n, n)`, then set the non-zero elements.
"""
struct System
    magnetic_field::MagneticField
    nuclei::Vector{Nucleus}
    chemical_shifts::Vector{Float64}
    couplings::Matrix{Frequency}

    function System(
        magnetic_field::MagneticField,
        nuclei::Vector{Nucleus},
        chemical_shifts::Vector{Float64},
        couplings::Matrix{<:Frequency},
    )
        n = length(nuclei)
        @assert length(chemical_shifts) == n
        @assert size(couplings) == (n, n)
        for i in 1:n
            @assert ustrip(couplings[i, i]) == 0
            for j in 1:i-1
                @assert ustrip(couplings[i, j]) == 0 || ustrip(couplings[j, i]) == 0
            end
        end
        new(magnetic_field, nuclei, chemical_shifts, couplings)
    end
end

end # module SpinSystem

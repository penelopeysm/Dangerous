module SpinSystem

using ..Nuclei: Nucleus

struct System
    magnetic_field::Float64  # TODO: specify units, right now it's MHz
    nuclei::Array{Nucleus}
    chemical_shifts::Array{Float64}
    couplings::Matrix{Float64}
end

end # module SpinSystem

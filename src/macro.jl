function add_args(expr)
    # rewrite f(a, b, c, ...) into 
    # ρ = f(ρ, sys, a, b, c, ...)
    #
    if expr.head != :call
        error("Expected a call expression")
    end
    # Try to get the function name from a lookup. Otherwise assume it's a
    # user-defined function and escape it.
    func = get(_MACRO_LOOKUP, expr.args[1], esc(expr.args[1]))
    return Expr(
        Symbol("="),
        :ρ,
        Expr(
            expr.head,           # This is just a symbol, don't escape it
            func,                # The function we're applying
            :ρ,                  # Don't escape this because it's defined in the macro
            esc(:sys),           # Escape this because it's passed in
            [esc(x) for x in expr.args[2:end]]...  # Escape the rest of the arguments
        )
    )
end

macro pulse_sequence(sys, sequence_expr)
    return Expr(:block,
        :(ρ = SpinSystemState(ρ_eq($(esc(sys))))),
        add_args(sequence_expr.args[2]),
        :(return collapse(ρ))
    )
end

"""
    pulse_instant(nuc, angle, phase)

NOTE: This method should only be used inside the @pulse_sequence macro.

Apply an instantaneous pulse to all `nuc` spins in the system, with a given
flip angle (in radians) and phase. The phase can be specified either as an
angle in radians or as the following symbols: `:x`, `:y`, `:_x`, or `:_y`,
which correspond respectively to 0, π/2, π, and 3π/2.
"""
function pulse_instant(nuc, angle, phase)
    error("This method should only be used inside the @pulse_sequence macro.")
    nuc, angle, phase  # To silence LSP
end

function _psmacro_pulse_instant(ρ, sys, nuc, angle, phase)
    return propagate(ρ, h_pulse(sys, nuc, angle * u"Hz", phase), 1u"s")
end

const _MACRO_LOOKUP = Dict(
    :pulse_instant => :_psmacro_pulse_instant
)

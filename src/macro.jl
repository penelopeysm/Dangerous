function add_args(expr)
    # rewrite f(a, b, c, ...) into 
    # ρ = f(ρ, sys, a, b, c, ...)
    #
    if expr.head != :call
        error("Expected a call expression")
    end
    return Expr(
        Symbol("="),
        :ρ,
        Expr(
            expr.head,           # This is just a symbol, don't escape it
            esc(expr.args[1]),   # Escape this because it's the function name
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

"""
    pulse_instant(ρ, sys, nuc, angle, phase)

Applies an instantaneous pulse (with parameters as above) to the spin system
`sys` with the state `ρ`.
"""
function pulse_instant(ρ, sys, nuc, angle, phase)
    return propagate(ρ, h_pulse(sys, nuc, angle * u"Hz", phase), 1u"s")
end

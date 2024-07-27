"""
Transforms f(a, b, c, ...) into f(ρ, sys, a, b, c, ...)
"""
function add_args_to_function(expr::Expr)
    if expr.head != :call
        error("Expected a call expression")
    end
    return Expr(
        expr.head,           # This is just a symbol, don't escape it
        esc(expr.args[1]),   # Escape this because it's the function name
        :ρ,                  # Don't escape this because it's defined in the macro
        esc(:sys),           # Escape this because it's passed in
        [esc(x) for x in expr.args[2:end]]...  # Escape the rest of the arguments
    )
end

function transform_expr(line_number_node::LineNumberNode)
    return [line_number_node]
end

function transform_expr(expr::Expr)
    if expr.head == :call
        # rewrite
        #     f(a, b, c, ...)
        # into 
        #     ρ = f(ρ, sys, a, b, c, ...)
        return [Expr(
            Symbol("="),
            :ρ,
            add_args_to_function(expr)
        )]
    elseif expr.head == Symbol("=")
        # rewrite
        #     var = f(a, b, c, ...)
        # into 
        #     global var = f(ρ, sys, a, b, c, ...)
        #     ρ = var.ρ
        expr1 = Expr(
            :global,
            Expr(
                Symbol("="),
                esc(expr.args[1]),    # Escape because we want the user to be able to get this value
                add_args_to_function(expr.args[2]),
            )
        )
        expr2 = :(ρ = $(esc(expr.args[1])).ρ)
        return [expr1, expr2]
    else
        error("Unsupported expression type: $(expr.head)")
    end
end

macro pulse_sequence(sys, sequence_expr)
    transformed_exprs = []
    for expr in sequence_expr.args
        push!(transformed_exprs, transform_expr(expr)...)
    end

    return Expr(:block,
        :(ρ = SpinSystemState(ρ_eq($(esc(sys))))),
        transformed_exprs...,
        # :(return collapse(ρ))  # Not necessary (yet)
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


"""
    delay(duration)

NOTE: This method should only be used inside the @pulse_sequence macro.

A delay to allow for free evolution.
"""
function delay(duration::Time)
    error("This method should only be used inside the @pulse_sequence macro.")
    duration  # To silence LSP
end


"""
    delay(ρ, sys, duration)

Propagate the state `ρ` of the spin system `sys` under free evolution for the
given duration.
"""
function delay(ρ, sys, duration::Time)
    return propagate(ρ, h_free(sys), duration)
end

using Dangerous
using Logging
using Unitful: @u_str

sys = System(
    14.1u"T",
    Dict(H1 => 4.7),
    [H1, H1],
    [1.5, 9],
    [0 30; 0 0]u"Hz"
)

# Either use the builtin functions

@pulse_sequence sys begin
    pulse_instant(H1, π / 2, :_y)
end

# You can define your own functions if you want, but there'll be LSP errors

function my_pulse_instant(ρ, sys, nuc, angle, phase)
    return Dangerous.propagate(ρ, Dangerous.h_pulse(sys, nuc, angle * u"Hz", phase), 1u"s")
end

@pulse_sequence sys begin
    my_pulse_instant(H1, π / 2, :_y)
end



@macroexpand @pulse_sequence sys begin
    pulse_instant(H1, π / 2, :_y)
    pulse_instant(H1, π, :x)
end

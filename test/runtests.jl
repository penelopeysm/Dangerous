import Dangerous as D
using Test: @test

maybe_Y = D.propagate(D.SingleSpin.Z, D.SingleSpin.h_pulse(π/2, π), 1)
@test isapprox(real.(maybe_Y), real.(D.SingleSpin.Y), atol=1e-15)
@test isapprox(imag.(maybe_Y), imag.(D.SingleSpin.Y), atol=1e-15)

.PHONY: docs
docs:
	julia --project=docs --color=yes docs/make.jl

.PHONY: test
test:
	julia --project=. -e "import Pkg; Pkg.test()"

.PHONY: docs test

docs:
	julia --project=docs --color=yes docs/make.jl

test:
	julia --project=. -e "import Pkg; Pkg.test()"

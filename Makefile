.PHONY: docs serve-docs test

docs:
	julia --project=docs --color=yes docs/make.jl

serve-docs:
	python -m http.server -d docs/build

test:
	julia --project=. -e "import Pkg; Pkg.test()"

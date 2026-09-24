JULIA:=julia

default: help

setup:
	${JULIA} --project=@runic --startup-file=no -e 'using Pkg; Pkg.add("Runic")'

format:
	${JULIA} --project=@runic --startup-file=no -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace .

testgroups:
	${JULIA} test/group_list.jl

test:
	${JULIA} --project -e 'using Pkg; Pkg.update(); Pkg.test(; test_args = ARGS)' -- $(ARGS)

docs:
	${JULIA} --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.update()'
	${JULIA} --project=docs docs/make.jl

all: setup format test docs

help:
	@echo "The following make commands are available:"
	@echo " - make setup: install the dependencies for make command"
	@echo " - make format: format codes with JuliaFormatter"
	@echo " - make testgroups: show the available test groups"
	@echo " - make test: run the tests"
	@echo " - make docs: instantiate and build the documentation"
	@echo " - make all: run every commands in the above order"

.PHONY: default setup format testgroups test docs all help

# Need to load AstrochemicalYields into Main to work with ParallelTestRunner
@eval Main using AstrochemicalYields
using Documenter: DocMeta, doctest

DocMeta.setdocmeta!(Main.AstrochemicalYields, :DocTestSetup, :(using AstrochemicalYields); recursive=true)
doctest(Main.AstrochemicalYields)

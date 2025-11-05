import AstrochemicalYields
using ParallelTestRunner: runtests, parse_args

# Can run single test files with, e.g., Pkg.test("AstrochemicalYields"; test_args=`--verbose doctests`)
# Can list available single tests with Pkg.test("AstrochemicalYields"; test_args=`--list`);
const init_code = quote
    using AstrochemicalYields
end
args = parse_args(Base.ARGS)
runtests(AstrochemicalYields, args; init_code)

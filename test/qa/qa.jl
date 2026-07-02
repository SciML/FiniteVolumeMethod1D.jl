using SciMLTesting, FiniteVolumeMethod1D, Test
using JET

run_qa(
    FiniteVolumeMethod1D;
    explicit_imports = true,
)

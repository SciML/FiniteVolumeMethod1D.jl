using SciMLTesting, FiniteVolumeMethod1D, Test
using JET

run_qa(
    FiniteVolumeMethod1D;
    explicit_imports = true,
    ei_kwargs = (;
        # Other packages' not-yet-public names; ignore until they go public.
        all_qualified_accesses_are_public = (;
            ignore = (
                :AutoSpecialize,  # SciMLBase.AutoSpecialize (specialization default)
                :init,            # CommonSolve.init (extended/forwarded here)
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (:solve,),   # CommonSolve.solve (re-exported)
        ),
    ),
)

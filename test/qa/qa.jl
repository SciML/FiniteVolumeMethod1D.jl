using SciMLTesting, FiniteVolumeMethod1D, Test
using JET

run_qa(
    FiniteVolumeMethod1D;
    api_docs_kwargs = (;
        rendered = true,
        rendered_ignore = (:solve,),
    ),
    explicit_imports = true,
)

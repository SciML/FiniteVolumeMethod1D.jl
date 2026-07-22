using SciMLTesting, FiniteVolumeMethod1D, Test
run_qa(FiniteVolumeMethod1D; reexports_allow = (:solve,))

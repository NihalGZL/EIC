# EIC

To be able compile the macros in this project we need to link to certain libraries first.

$  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/nihal/pythia8310/lib 

$  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/nihal/pythia8310/hepmc3-install/lib

$  export PYTHIA8DATA=/home/nihal/pythia8310/share/Pythia8/xmldoc

 Only difference between EEC.C and generate_hepmc.C is generate_hepmc.C is saving events in HEPMC format.


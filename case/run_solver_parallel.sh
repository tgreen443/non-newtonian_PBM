 #!/bin/bash

nprocs=32
foamDictionary system/decomposeParDict -entry numberOfSubdomains -set $nprocs

rm -rf 0

cp -r 0_org 0

cp 0/alpha.air.org 0/alpha.air
cp 0/alpha.water.org 0/alpha.water

setFields | tee log.setFields

decomposePar | tee log.decomposePar

mpirun -np $nprocs renumberMesh -overwrite -parallel  | tee log.renumberMesh

mpirun -np $nprocs multiphaseEulerFoam -parallel | tee log.solver

reconstructPar | tee log.reconstructPar




#When running with two-phases the optional referencePhase entry in phaseProperties can be used to specify
#which phase fraction should not be solved, providing compatibility with reactingtTwoPhaseEulerFoam.


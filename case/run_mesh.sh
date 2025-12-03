#!/bin/bash

rm -rf 0

fluent3DMeshToFoam 1.msh | tee log.fluent3DMeshToFoam

checkMesh | tee log.checkMesh

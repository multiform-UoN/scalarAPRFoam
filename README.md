# scalarAPRFoam

A scalar reactive transport solver for OpenFOAM 8, with a surface reaction boundary condition.

This solver was developed by Matteo Icardi and Diego Fida.

## Description

`scalarAPRFoam` is a transient solver for multiple species transport with arbitrary volumetric and surface reactions. The coupling between species is implemented using a segregated algorithm.

The solver is based on the `pimpleFoam` solver and includes:
-   Transport equations for multiple species.
-   Volumetric reactions.
-   A surface reaction boundary condition.

## Notes

-   This solver is developed for OpenFOAM 8.
-   A scientific paper describing the solver and its application is in preparation. The authors of the paper include Matteo Icardi, Diego Fida, and others.
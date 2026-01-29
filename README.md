# scalarAPRFoam

A scalar reactive transport solver for OpenFOAM available for OpenFOAM-8 and OpenFOAM-13, with a surface reaction boundary condition now available in OpenFOAM-8.

This solver was developed by Matteo Icardi and Diego Fida.

## Description

`scalarAPRFoam` is a transient solver for multiple species transport with arbitrary volumetric and surface reactions. The coupling between species is implemented using a segregated algorithm.

The solver is based on the `pimpleFoam` solver and includes:
-   Transport equations for multiple species.
-   Volumetric reactions.
-   A surface reaction boundary condition (now available only on OpenFOAM-8).

The transport equation is :

$\frac{\partial\left(\varepsilon_{\text{cell}}\ c_i\right)}{\partial t} + \mathbf{u} \cdot \nabla c_i = \varepsilon_{\text{cell}}D_{i} \nabla^2 c_i + \varepsilon_{\text{cell}}R_i\ ,$

where the $\varepsilon_{\text{cell}}$ is the cell porosity defined by cellsZone, $c_i$ the concentratio, $\mathbf{u}$ the velocity, $D_i$ the diffusion coefficient defined as an armonic average of coefficient in the fluid and in the solid, finaly, $R_i$ is the reaction reate. 

The kinetic equation available for the reaction rate are rappresented by:

$R_{i} =\frac{k_i\prod^n_{j=1}c_j^{a_{j,i}}}{\left(1+\sum^n_{j=1}\left(K_{j,i}c_j\right)^{b_{j,i}}\right)^{c_i}}$

and is implemented by splitting this term in an implicit and explicit contribution:

$R_i=  R_i(\textbf{c}^0) + \sum_{j=1}^n \frac{\partial R_i(\textbf{c}^0)}{\partial c_j} (c_j-c_j^0)$

$K=  \frac{\partial R_i(\textbf{c}^0)}{\partial c_i}$

$F =    R_i(\mathbf{c}^0) - \frac{\partial R_i(\textbf{c}^0)}{\partial c_i}  \ c_i^0 $ 

## Notes
-   A scientific paper describing the solver and its application is in preparation. The authors of the paper include Matteo Icardi, Diego Fida, and others.
-   The tutorial present is developed in OpenFoam-8 
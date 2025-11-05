/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 \\    /   O peration     | Website:  https://openfoam.org
  \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
   \\/     M anipulation  |
-------------------------------------------------------------------------------
License
This file is part of OpenFOAM.

OpenFOAM is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
reactiveTransportFoam

Description
Transient solver for multiple species transport with
arbitrary volumetric reaction
Coupling is implemented using a segregated algorithm.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define specie from dictionary
typedef dictionary specie;

int main(int argc, char *argv[])
{
  #include "postProcess.H"
  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"

  //Add pimple coupling controls
  pimpleControl pimple(mesh);

  #include "createFields.H"
  #include "CourantNo.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


  //*** search the reactive zone ***//

  label zoneIndex = -1;                           //- index of the reactive and solid zone 

  if (reactiveZoneName=="all")                    //- if the reactive zone is all the domain  flagReactive is set to 1
  {
    flagReactive = 1.0;
    Info<< "All the domain is reactive" << endl;
  }
  else                                            //- if the reactive zone is a specific zone, flagReactive is set to 1 only in the reactive zone
  {
    zoneIndex = mesh.cellZones().findZoneID(reactiveZoneName);
    if (zoneIndex==-1)                            //- if the zone does not exist, a warning is printed
    {
      Foam::Warning();
      Info << "The zone " << reactiveZoneName << " does not exist" << endl;
    }
    else
    {
      const cellZone& targetZone = mesh.cellZones()[zoneIndex];
      forAll(targetZone, cellI)                   //- loop all the cell to define the field flagReactive
      {
          flagReactive[targetZone[cellI]] = 1.0;
      }
    }
  }
  
  forAll(reaction, re)                          //- loop over all the reactions to define the different flagReactive field for each reaction
  {
    flagReactivePtrL.set(
      re,
      new volScalarField
      (
          IOobject
          (
              "flagReactive_"+reactionEntries[re].keyword(),
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedScalar
          (
              "zero",
              dimensionSet(0, 0, 0, 0, 0),
              0.0
          )
      )
    );  

    if (reaction[re].found("zoneName"))               //- if a specific reactive zone is defined for the reaction it overwritten the global one, else the global one is used
    {
      Info<< "Reaction " << reactionEntries[re].keyword() << " set a specific reactive zone" << endl; 
      word reactiveZoneName(reaction[re].lookup("zoneName"));
      zoneIndex = -1;
      volScalarField& specificFlag = flagReactivePtrL[re];

      if (reactiveZoneName=="all")                    //- if the reactive zone is all the domain  flagReactive is set to 1
      {
        specificFlag = 1.0;
      }
      else                                            //- if the reactive zone is a specific zone, flagReactive is set to 1 only in the reactive zone 
      {
        zoneIndex = mesh.cellZones().findZoneID(reactiveZoneName);
        if (zoneIndex==-1)                            //- if the zone does not exist, a warning is printed
        {
          Foam::Warning();
          Info << "The zone " << reactiveZoneName << " selected for " << reactionEntries[re].keyword()<< " does not exist" << endl;
        }
        else
        {
          const cellZone& targetZone = mesh.cellZones()[zoneIndex];
          forAll(targetZone, cellI)                   //- loop all the cell to define the field flagReactive
          {
            specificFlag[targetZone[cellI]] = 1.0;
          }
        }
      }
    }
    else
    {
      flagReactivePtrL[re] = flagReactive;
    }
  }
  


 
  //*** create effective diffusivity surface field ***//

  forAll(species, sp)                               //- loop over all the species                    
  {
    dimensionedScalar D(species[sp].lookup("D"));   //- Diffusivity in the fluid  
    dimensionedScalar Ds(species[sp].lookup("D"));  //- Diffusivity in the solid initializated with the fluid value
    
    if (zoneIndex==-1)                              //- If the zone does not exist, the solid diffusivity is the same as the fluid
    {
      Info<< "Solid zone does not selected" << endl;
    }
    else                                          //- If the zone exist, the solid diffusivity is defined in the dictionary
    {
      Ds = species[sp].lookup("Ds");
    }


    //- Create effective diffusivity SURFACE field as a armonic mean of VOL diffusivity, to take into account the different diffusivity in the solid and in the fluid
    DsupPtrL.set
    (
      sp,
      new surfaceScalarField
      (
        IOobject
        (
          "Dsup_" + specEntries[sp].keyword(),
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        // harmonic mean on faces
        1.0 / reverseLinear<scalar>(mesh).interpolate
        (
            ((scalar(1) - flagReactive)/D + flagReactive/Ds)
        )
      )
    );
  }



  
  Info<< "\nStarting time loop\n" << endl;

  while (pimple.loop(runTime))
  {
    Info<< "Time = " << runTime.timeName() << nl << endl;
   
    while (pimple.loop())                         //- Coupling loop
    {

      #include "updateReactionRates.H"

      forAll(species,sp)                          //- loop over all the species
      {

        //Gather pointers to specific species sp
        volScalarField&   C(cPtrL[sp]);
        volScalarField&   Rk(RkPtrL[sp]);
        volScalarField&   R(RPtrL[sp]);
        surfaceScalarField& Dsup(DsupPtrL[sp]);


        while (pimple.correctNonOrthogonal())     //- Non-orthogonal correction loop
        {

          fvScalarMatrix CEqn
          (
            fvm::ddt(C)
            + fvm::div(phi,C,"div(phi,C)")
            - fvm::laplacian(Dsup,C,"laplacian(D,C)")
            ==
            - fvm::SuSp(-Rk,C)                      //- Add the implicit part to the matrix diagonal if negative
            + R                                     //- Add total reaction rate as explicit source term
            - Rk*C                                  //- Subtract the implicit part as explicit source term  to avoid double counting
            + fvOptions(C)
          );

          CEqn.relax();
          fvOptions.constrain(CEqn);
          CEqn.solve();
          fvOptions.correct(C);
        }

      } //end loop species

    } // end pimple loop

    runTime.write();

  }

  Info<< "End\n" << endl;
  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;


  return 0;
}


// ************************************************************************* //

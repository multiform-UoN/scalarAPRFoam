
Info<< "Reading field U\n" << endl;
volVectorField U                        //- Create Velocity field U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "velocity",
        dimensionSet(0,1,-1,0,0,0,0),
        vector(0,0,0)
    )
);


Info<< "Reading reaction dictionary\n" << endl;
IOdictionary reactionDict           //- Create reactionDict dictionary to define the reactions
(
    IOobject
        (
            "reactionDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
);

PtrList<dictionary> reaction;
PtrList<entry> reactionEntries
(
    reactionDict.lookup("reactions")
);
scalar porosity = reactionDict.lookup<scalar>("porosity");
reaction.setSize(reactionEntries.size());
PtrList<volScalarField> rPtrL(reaction.size());
PtrList<volScalarField> flagReactivePtrL(reaction.size());

forAll(reaction, re)
{
    reaction.set(re, new dictionary(reactionEntries[re].dict()));

    rPtrL.set(
        re,
        new volScalarField
        (
            IOobject
            (
                "r_"+reactionEntries[re].keyword(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, -3, -1, 0, 1, 0, 0),
                0.0
            )
        )
    );
}


Info<< "Reading reactiveZone name\n" << endl;
word reactiveZoneName(reactionDict.lookup("zoneName"));             //- Name of the reactive zone

volScalarField flagReactive                                         //- Flag field to identifie reactive cells=1 and inert cells=0
(
    IOobject
    (
        "flagReactive",
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
);


Info<< "Reading species dictionary\n" << endl;
IOdictionary speciesDict                                            //- dictionary to define the species
(
    IOobject
        (
            "speciesDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
);

//- Read chemical species
PtrList<dictionary> species;
PtrList<entry> specEntries(speciesDict.lookup("species"));
species.setSize(specEntries.size());

Info<< "Reading/calculating species properties\n" << endl;


PtrList<volScalarField> cPtrL(species.size());                  //- Create concentration fields


PtrList<volScalarField> RkPtrL(species.size());                 //- Create implicit reaction rates fields
PtrList<volScalarField> RPtrL(species.size());                  //- Create total reaction rates fields
PtrList<surfaceScalarField> DsupPtrL(species.size());           //- Create effective diffusivity fields

forAll(species,sp)                                              //- loop over all species
{
    
    species.set(sp,new dictionary(specEntries[sp].dict()));     //- Create a dictionary for each species

    cPtrL.set(                                                  //- Create concentration fields for i-th species
        sp,
        new volScalarField
        (
            IOobject
            (
                specEntries[sp].keyword(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );


    RkPtrL.set(                                                 //- Create implicit reaction rates for i-th species
        sp,
        new volScalarField
        (
            IOobject
            (
                "Rk_"+specEntries[sp].keyword(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, 0, -1, 0, 0, 0, 0),
                0.0
            )
        )
    );

    RPtrL.set(                                                   //- Create total reaction rates for i-th species
        sp,
        new volScalarField
        (
            IOobject
            (
                "R_"+specEntries[sp].keyword(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, -3, -1, 0, 1, 0, 0),
                0.0
            )
        )
    );

}


dimensionedScalar SMALL_C                                         //- Smallest machine number with concentration dimensions   
(
    "SMALL_C", 
    dimensionSet(0, -3, 0, 0, 1, 0, 0),
    SMALL
);


#include "createFvOptions.H"
#include "createPhi.H"

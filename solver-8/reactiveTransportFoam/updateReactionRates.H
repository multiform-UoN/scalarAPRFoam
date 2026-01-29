//              k*prod(C_i^exponents_i)
// r = ------------------------------------------------
//     {1+sum[(Ks_i*C_i)^exponentsS_i]}^exponentsS_total





forAll(species, j)                      //- Clean the reaction rates
{
    volScalarField& Rspec = RPtrL[j];
    Rspec = dimensionedScalar("zero", Rspec.dimensions(), 0.0);

    volScalarField& Rkspec = RkPtrL[j];
    Rkspec = dimensionedScalar("zero", Rkspec.dimensions(), 0.0);
}


forAll(reaction, re)                    //- loop over all the reactions
{

    dimensionedScalar k (reaction[re].lookup("k"));                            //- Reaction rate constant
    const Foam::wordList& reactiveSpecies = reaction[re].lookup("species");    //- Species involved in the reaction
    const scalarField& sto = reaction[re].lookup("stoichiometry");             //- Stoichiometry of the reaction
    const scalarField& n = reaction[re].lookup("exponents");                   //- Exponents of each species in kinetic law
    const scalarField& kS = reaction[re].lookup("kS");                         //- Adsorption rate constant
    const scalarField& ns = reaction[re].lookup("exponentsS");                 //- Exponents of each species in adsorption law
    const label ns_total = reaction[re].lookupOrDefault<label>("exponentsS_total", 1);     //- Total exponent of the adsorption law (adimensional scalar)
    List<label> i(reactiveSpecies.size(), -1);                                 //- Index of the reactive species
    volScalarField& thisReactionFlag = flagReactivePtrL[re];                   //- Reactive flag field


    volScalarField& r = rPtrL[re];                                             //- This reaction rate
    forAll (r, cellI)
    {
        r[cellI] = k.value() * thisReactionFlag[cellI];                         //- Reaction rate initialization to k inside the reactive zone and zero outside
    }

    volScalarField reactionDenominator                                          //- Reaction rate denominator initialized to one
    (
        IOobject
        (
            "reactionDenominator",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)
    );

    forAll(reactiveSpecies, j)          //- loop over all the reactive species
    {
        bool found = false;

        forAll(species, sp)             //- loop over all species to correlate the position of the species in the species dictionary with the position in the reactive species list
        {

            if (reactiveSpecies[j] == specEntries[sp].keyword())
            {
                i[j] = sp;
                found = true;
                break;
            }
        }

        if (!found)                     //- print a warning if the reaction dictionary contains a species that does not exist
        {
            Foam::Warning();
            Info<< "The species " << reactiveSpecies[j] << " does not exist" << endl;
        }

    }

    forAll(reactiveSpecies, j)            //- loop over all the reactive species to calculate the reaction rate
    {
        if (i[j] >= 0)                    //- leaves away -1, the species that do not exist
        {
            volScalarField& C = cPtrL[i[j]];        //- Concentration field

            if (n[j] > 0)
            {
                r *= Foam::pow(C*Foam::pos(C),n[j])/dimensionedScalar("one", dimensionSet(0, -3*n[j], 0, 0, n[j], 0, 0), 1.0);                             //- Update the reaction rate multiplying by the "n" power of the current species
            }
            
            else if(n[j] < 0)
            {
                // Add a small constant to avoid division by zero when C=0 and n<0
                r *= Foam::pow((C*Foam::pos(C))+SMALL_C,n[j])*Foam::pos(C)/dimensionedScalar("one", dimensionSet(0, -3*n[j], 0, 0, n[j], 0, 0), 1.0);      //- Update the reaction rate, but avoid the ratio 1/0 when n<1 and C=0
            }

            reactionDenominator += Foam::pow(kS[j]*C*Foam::pos(C), ns[j])/dimensionedScalar("one", dimensionSet(0, -3*ns[j], 0, 0, ns[j], 0, 0), 1.0);         //- Update the reaction rate denominator adding the "nS_j" power of the current species adsorption
        }
    }
        r /= pow(reactionDenominator, ns_total);                                                    //- Update the reaction rate by the denominator

    forAll(reactiveSpecies, j)
    {
        if (i[j] >= 0)
        {
            volScalarField& C = cPtrL[i[j]];                                     //- Concentration field of the particular species

            volScalarField& Rspec = RPtrL[i[j]];                                 //- Reaction rate field of the particular species
            Rspec += sto[j] * r;                                                 //- Update the reaction rate field of the particular species
            
            volScalarField& Rkspec = RkPtrL[i[j]];                               //- Implicit reaction rate field of the particular species
            Rkspec += sto[j] * r * ( n[j] / (C+SMALL_C)*Foam::pos(C) - ns_total*ns[j]*Foam::pow(kS[j], ns[j])*Foam::pow(C*Foam::pos(C)+SMALL_C , ns[j]-1)*Foam::pos(C)/reactionDenominator/dimensionedScalar("one", dimensionSet(0, -3*ns[j], 0, 0, ns[j], 0, 0), 1.0)); //- Update the implicit reaction rate field of the particular species
        }
    }
}

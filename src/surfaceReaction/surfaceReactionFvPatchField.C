/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 The OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                      Copyright (C) 2025 Matteo Icardi, Diego Fida
-------------------------------------------------------------------------------
License
    This file is part of scalarAPRFoam.

    scalarAPRFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    scalarAPRFoam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with scalarAPRFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    This boundary condition implements a Robin boundary condition for a surface
    reaction.

    The reaction rate is given by a general kinetic law, which can include an
    adsorption term. The parameters of the reaction (rate constants, stoichiometric
    coefficients, exponents, etc.) are read from the `surfaceReactionDict` dictionary.

\*---------------------------------------------------------------------------*/

#include "surfaceReactionFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceReactionFvPatchField::surfaceReactionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(p, iF),
    dict_(new dictionary()),
    dictName_("surfaceReactionDict"),
    RobinKeff_(p.size(), 0.0),
    RobinFeff_(p.size(), 0.0)
{}


Foam::surfaceReactionFvPatchField::surfaceReactionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    RobinFvPatchScalarField(p, iF, dict),
    dict_(dict),
    dictName_(dict.lookupOrDefault<word>("reactionDict", "surfaceReactionDict")),
    // dict_(dict.subDict("reactions")),
    RobinKeff_(p.size(), 0.0),
    RobinFeff_(p.size(), 0.0)
{
    IOdictionary reactionDict
                (
                    IOobject
                        (
                            dictName_,
                            internalField().mesh().time().constant(),
                            internalField().mesh(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                );
    dict_ = reactionDict;
}


Foam::surfaceReactionFvPatchField::surfaceReactionFvPatchField
(
    const surfaceReactionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    RobinFvPatchScalarField(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    dictName_(ptf.dictName_),
    RobinKeff_(mapper(ptf.RobinKeff_)),
    RobinFeff_(mapper(ptf.RobinFeff_))
{}


Foam::surfaceReactionFvPatchField::surfaceReactionFvPatchField
(
    const surfaceReactionFvPatchField& ptf
)
:
    RobinFvPatchScalarField(ptf),
    dict_(ptf.dict_),
    dictName_(ptf.dictName_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_)
{}


Foam::surfaceReactionFvPatchField::surfaceReactionFvPatchField
(
    const surfaceReactionFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(ptf, iF),
    dict_(ptf.dict_),
    dictName_(ptf.dictName_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::surfaceReactionFvPatchField::autoMap
(
    const fvPatchFieldMapper&   m
)
{
    RobinFvPatchScalarField::autoMap(m);
    m(RobinKeff_,RobinKeff_);
    m(RobinFeff_,RobinFeff_);
}

void Foam::surfaceReactionFvPatchField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    RobinFvPatchScalarField::rmap(ptf,addr);

    const surfaceReactionFvPatchField& mptf =
        refCast<const surfaceReactionFvPatchField>(ptf);

    RobinKeff_.rmap(mptf.RobinKeff_,addr);
    RobinFeff_.rmap(mptf.RobinFeff_,addr);
}

void Foam::surfaceReactionFvPatchField::write(Ostream& os) const
{
    RobinFvPatchScalarField::write(os);
    // writeEntry(os, "reactions", dict_);
    writeEntry(os, "reactionDict", dictName_);
    writeEntry(os, "RobinKeff", RobinKeff_);
    writeEntry(os, "RobinFeff", RobinFeff_);
}

void Foam::surfaceReactionFvPatchField::evaluate
(
    const Pstream::commsTypes commsType
)
{

    updateCoeffs();

    //- Evaluate Robin boundary condition
    RobinFvPatchScalarField::evaluate();

}


void Foam::surfaceReactionFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
      return;
    }

    // Internal loop coefficients initialization
    const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
    const scalarField& RobinF = RobinFvPatchScalarField::RobinF();

    // Effective reaction rate coefficients initialization
    RobinFeff_ = RobinF;
    RobinKeff_ = RobinK;

    // Iterate over all entries in the parent dictionary (each iteration is a reaction)
    forAll(dict_.toc(), i)
    {

        const dictionary& entry = dict_.subDict(dict_.toc()[i]);

        const scalar k(entry.lookup<scalar>("k"));                  // Reaction rate constant

        const Foam::wordList& species = entry.lookup("species");    // Species involved in the reaction

        const scalarField& sto = entry.lookup("stoichiometry");     // Stoichiometry of the reaction

        const scalarField& a = entry.lookup("exponents");           // Exponents of each species in kinetic law

        const scalarField& kS = entry.lookup("kS");                 // Adsorption rate constant

        const scalarField& b = entry.lookup("exponentsS");          // Exponents of each species in adsorption law

        const label exponentsS_total = entry.lookupOrDefault<label>("exponentsS_total", 1); // Total exponent of the adsorption law (adimensional scalar)

        label j0(-1);                                                // Initialization index of the current species

        // Calculate the reaction rate based on the kinetic law
        //
        // r = k * product(C_i^a_i) / (1 + sum(kS_j*C_j^b_j))^exponentsS_total
        //
        scalarField reaction(this->size(), k);                      // Reaction rate initialization to k 
        
        forAll(species, j)                                          // Loop over all species to calculate the numerator of the reaction rate
        {
            if (species[j] == this->internalField().name())         // If the current species is the field solved by the BC
            {
                j0 = j;                                             // Save the index of the current species
                reaction *= Foam::pow(*this*Foam::pos(*this),a[j]); // Update the reaction rate multiplying by the "a" power of the current species
            }
            else
            {
                const scalarField& Yj = patch().lookupPatchField<volScalarField, scalar>(species[j]); // Lookup the internal loop species
                reaction *= Foam::pow(Yj*Foam::pos(Yj),a[j]);       // Update the reaction rate multiplying by the "a_j" power of the current species
            }
        }                                                           // End loop over all species for numerator r=k*prod(Y_i^exponents_i)

        
        scalarField reactionDenominator(this->size(), 1.0);         // Reaction rate denominator initialization to 1 (1+sum^n_{j=1}[kS_{j}Y_j]^{exponentsS_j})**{exponentsS_total}

        forAll(species, j)                                          // Loop over all species to calculate the denominator of the reaction rate
        {
            if (b[j]>0)
            {
                if (species[j] == this->internalField().name())         // If the current species is the field solved by the BC
                {
                    reactionDenominator += Foam::pow(*this*Foam::pos(*this)*kS[j],b[j]);  // Update the reaction rate denominator adding the "b_j" power of the current species
                }
                else
                {
                    const scalarField& Yj = patch().lookupPatchField<volScalarField, scalar>(species[j]); // Lookup the internal loop species
                    reactionDenominator += Foam::pow(Yj*Foam::pos(Yj)*kS[j],b[j]);        // Update the reaction rate denominator adding the "b_j" power of the current species
                }
            }
        }

        reaction /= Foam::pow(reactionDenominator,exponentsS_total);                            // Update the reaction rate by the denominator

        // Linearize the reaction rate to obtain the coefficients of the Robin boundary condition
        //
        // R = R(C_0) + (C - C_0)*dR/dC|_{C_0} = F + K*(C - C_0)
        //
        // where:
        // F = R(C_0) - C_0*dR/dC|_{C_0}
        // K = dR/dC|_{C_0}
        //
        // This is done to improve the stability of the numerical solution.
        if (j0 >= 0)                                                // If the species was found do: else the species doesn't involve in the reaction 
        {
            // RobinKeff_ update summing the implicit part of this reaction (the implicit part is equal to the first concentration derivative)
            RobinKeff_ += sto[j0]*reaction * ( a[j0]/(*this + small) -exponentsS_total*b[j0]*Foam::pow(kS[j0],b[j0])*Foam::pow((*this*Foam::pos(*this)+small),b[j0]-1)/reactionDenominator);
            RobinFeff_ += sto[j0]*reaction;                         // RobinFeff_ update summing the reaction rate (implicit + explicit)
        }
    
    }

    RobinFeff_ -= RobinKeff_ *(*this);                              // RobinFeff_(explicit term) is updated by the subtraction of the explicit term 

    RobinFvPatchScalarField::updateCoeffs();                        // Evaluate Robin boundary condition


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        surfaceReactionFvPatchField
    );
}

// ************************************************************************* //

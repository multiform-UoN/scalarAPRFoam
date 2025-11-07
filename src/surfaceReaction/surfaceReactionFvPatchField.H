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

Class
    Foam::surfaceReactionFvPatchField

Description
    This boundary condition specialise the Robin boundary condition for
    an arbitrary surface reaction read from the dictionary surfaceReactionDict. This boundary is developed to a general kinetic law:
    
    $R_i=\sum_m \nu_{i,m} r_{i,m}$
    $r_{i,m} =\frac{k\prod^n_{j=1}Y_j^{exponents_j}}{(1+\sum^n_{j=1}kS_{j}Y_j^{exponentsS_j})^{exponentsS_total}}$

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceReaction;
        value           uniform 0;
    }
    \endverbatim


SourceFiles
    surfaceReactionFvPatchField.C

Developers:
    Matteo Icardi, Diego Fida (2024)
\*---------------------------------------------------------------------------*/


#ifndef surfaceReactionFvPatchField_H
#define surfaceReactionFvPatchField_H

#include "RobinFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class surfaceReactionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class surfaceReactionFvPatchField
:
    public RobinFvPatchScalarField
{
private:

    //- Dictionary
    dictionary dict_;
    word dictName_;

    //- Effective RobinK_
    scalarField RobinKeff_;

    //- Effective RobinF_
    scalarField RobinFeff_;

public:

    //- Runtime type information
    TypeName("surfaceReaction");


    // Constructors

        //- Construct from patch and internal field
        surfaceReactionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceReactionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given surfaceReactionFvPatchField
        //  onto a new patch
        surfaceReactionFvPatchField
        (
            const surfaceReactionFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceReactionFvPatchField
        (
            const surfaceReactionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceReactionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceReactionFvPatchField
        (
            const surfaceReactionFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceReactionFvPatchField(*this, iF)
            );
        }


        // Mapping functions

             //- Map (and resize as needed) from self given a mapping object
             virtual void autoMap
             (
                 const fvPatchFieldMapper&
             );

             //- Reverse map the given fvPatchField onto this fvPatchField
             virtual void rmap
             (
                 const fvPatchField<scalar>&,
                 const labelList&
             );


        // Member functions

            virtual scalarField& RobinK()
            {
                return RobinKeff_;
            }

            virtual const scalarField& RobinK() const
            {
                return RobinKeff_;
            }


            virtual scalarField& RobinF()
            {
                return RobinFeff_;
            }

            virtual const scalarField& RobinF() const
            {
                return RobinFeff_;
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

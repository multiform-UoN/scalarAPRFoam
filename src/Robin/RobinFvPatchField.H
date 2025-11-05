/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) The OpenFOAM Foundation, Ltd.
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
  
  Modifications and additional contributions:
  
  MultiForm Group, University of Nottingham
  Copyright (C) Matteo Icardi and collaborators
  Contributions to this file are licensed under the same GPLv3 terms as OpenFOAM.

  This work is based on OpenFOAM, with substantial portions copied, modified, or 
  extended under the GNU General Public License (GPLv3).
  Please refer to the original OpenFOAM copyright notice for the base framework, 
  and to MultiForm Group for new additions or modifications. For the full terms of 
  this license, see <http://www.gnu.org/licenses/>.



Class
    Foam::RobinFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a base class for 'Robin' type boundary
    conditions
        \f[
            D \nabla_\perp x  = K x + F
        \f]

    where
    \vartable
        x_p   | patch values
        x_c   | patch internal cell values
        K     | reaction term
        \Delta| inverse distance from face centre to internal cell centre
        F     | Forcing term
        D     | Diffusion
    \endvartable


Usage
    \table
        Property     | Description             | Required    | Default value
        RobinD       | Diffusion               | yes         |
        RobinF       | Source Term             | yes         |
        RobinR       | Reaction                | yes         |
    \endtable

SourceFiles
    RobinFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef RobinFvPatchField_H
#define RobinFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RobinFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RobinFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        //- Diffusion field
        scalarField RobinD_;

        //- Reaction field
        scalarField RobinK_;

        //- Source term
        Field<Type> RobinF_;

public:

    //- Runtime type information
    TypeName("Robin");


    // Constructors

        //- Construct from patch and internal field
        RobinFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RobinFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given RobinFvPatchField onto a new patch
        RobinFvPatchField
        (
            const RobinFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy
        RobinFvPatchField
        (
            const RobinFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new RobinFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        RobinFvPatchField
        (
            const RobinFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new RobinFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            virtual bool fixesValue() const
            {
                // True if RobinK is not zero
                forAll(RobinK_, i)
                {
                    if (abs(RobinK_[i]) > Foam::small)
                    {
                        return true;
                    }
                }
                return false;
            }

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


        // Return defining fields

            virtual scalarField& RobinD()
            {
                return RobinD_;
            }

            virtual const scalarField& RobinD() const
            {
                return RobinD_;
            }

            virtual scalarField& RobinK()
            {
                return RobinK_;
            }

            virtual const scalarField& RobinK() const
            {
                return RobinK_;
            }

            virtual Field<Type>& RobinF()
            {
                return RobinF_;
            }

            virtual const Field<Type>& RobinF() const
            {
                return RobinF_;
            }

        // Setting fields

            virtual void RobinD(scalarField& f)
            {
                RobinD_ = f;
            }

            virtual void RobinK(scalarField& f)
            {
                RobinK_ = f;
            }

            virtual void RobinF(Field<Type>& f)
            {
                RobinF_ = f;
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
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;

            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RobinFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

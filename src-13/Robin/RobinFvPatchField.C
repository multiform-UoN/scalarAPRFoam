/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 The OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                      Copyright (C) 2025 Matteo Icardi
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

\*---------------------------------------------------------------------------*/

#include "RobinFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    RobinD_(p.size(),1.0),
    RobinK_(p.size(),0.0),
    RobinF_(p.size(),pTraits<Type>::zero)
{
}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    RobinD_(p.size(),1.0),
    RobinK_(p.size(),0.0),
    RobinF_(p.size(),pTraits<Type>::zero)
{
  if (dict.found("RobinD"))
  {
    RobinD_ = scalarField("RobinD",dict,p.size());
  }
  else
  {
    RobinD_ = scalarField(p.size(),1.0);
  }
  if (dict.found("RobinK"))
  {
    RobinK_ = scalarField("RobinK",dict,p.size());
  }
  else
  {
    RobinK_ = scalarField(p.size(),0.0);
  }
  if (dict.found("RobinF"))
  {
    RobinF_ = Field<Type>("RobinF",dict,p.size());
  }
  else
  {
    RobinF_ = Field<Type>(p.size(),pTraits<Type>::zero);
  }
//    updateCoeffs();
//    evaluate();
}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
    RobinD_(mapper(ptf.RobinD_)),
    RobinK_(mapper(ptf.RobinK_)),
    RobinF_(mapper(ptf.RobinF_))
{
}

// ...existing code...

template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::RobinFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& m
)
{
    fvPatchField<Type>::map(ptf, m);
    
    const RobinFvPatchField<Type>& mptf =
        refCast<const RobinFvPatchField<Type>>(ptf);

    m(RobinD_, mptf.RobinD_);
    m(RobinK_, mptf.RobinK_);
    m(RobinF_, mptf.RobinF_);
}


template<class Type>
void Foam::RobinFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::reset(ptf);
    
    const RobinFvPatchField<Type>& mptf =
        refCast<const RobinFvPatchField<Type>>(ptf);

    RobinD_ = mptf.RobinD_;
    RobinK_ = mptf.RobinK_;
    RobinF_ = mptf.RobinF_;
}

template<class Type>
void Foam::RobinFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        (
          this->patchInternalField() * RobinD() * this->patch().deltaCoeffs()
          + RobinF()
        )
        /
        (
           RobinD() * this->patch().deltaCoeffs()
          - RobinK()
        )
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::snGrad() const
{
    return
        (
          RobinK() * this->patchInternalField()
          + RobinF()
        ) * this->patch().deltaCoeffs()
        /
        (
          RobinD() * this->patch().deltaCoeffs()
         - RobinK()
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one) *
    (
        RobinD() * this->patch().deltaCoeffs() )
        /
        (
            RobinD()*this->patch().deltaCoeffs() - RobinK()
        );

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         (
           RobinF()         /
           (
            RobinD() * this->patch().deltaCoeffs()
           - RobinK()
           )
         );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::gradientInternalCoeffs() const
{
    return Type(pTraits<Type>::one) *
    (
      RobinK()*this->patch().deltaCoeffs()
      /
      (
        RobinD() * this->patch().deltaCoeffs()
       - RobinK()
      )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
       (
         RobinF()* this->patch().deltaCoeffs()
         /
         (
          RobinD() * this->patch().deltaCoeffs()
         - RobinK()
         )
       );
}


template<class Type>
void Foam::RobinFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os,"RobinD",RobinD_);//RobinD_.writeEntry("RobinD", os);
    writeEntry(os,"RobinK",RobinK_);//RobinK_.writeEntry("RobinK", os);
    writeEntry(os,"RobinF",RobinF_);//RobinF_.writeEntry("RobinF", os);
    writeEntry(os,"value",*this);//this->writeEntry("value", os);
}



template<class Type>
void Foam::RobinFvPatchField<Type>::updateCoeffs()
 {

     if (this->updated())
     {
         return;
     }

     fvPatchField<Type>::updateCoeffs();
 }

// ************************************************************************* //

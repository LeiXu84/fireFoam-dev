/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "turbulentFluidThermoModel.H"

#include "mappedPatchBase.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::
nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("none"),
    massFluxFraction_(p.size(),0.)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::
nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    massFluxFraction_("massFluxFraction", dict, p.size())
{
// Ankur.. Is there a reason to NOT call mixedFvPatch constructor with (p, IF, dict) arguments here..??
    refValue() = 1.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }

    //if (dict.found("massFluxFraction"))
    //{
    //    massFluxFraction_ = Field<scalar>("massFluxFraction", dict, p.size());
    //}
    //else
    //{
    //    massFluxFraction_ = 0.;
    //}
}

Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::
nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    massFluxFraction_(ptf.massFluxFraction_,mapper)
{}


Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::
nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}

Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::
nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
// Ankur.. why don't we call autoMap function for mixedFvPatchField here.. ??
    scalarField::autoMap(m);
    massFluxFraction_.autoMap(m); 
}


void Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField& mptf = 
         refCast<const nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(ptf);
    
    massFluxFraction_.rmap(mptf.massFluxFraction_, addr);
}


void Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const label patchI = patch().index();

    const LESModel<EddyDiffusivity<compressible::turbulenceModel>>& turbModel =
        db().lookupObject
        <
            LESModel<EddyDiffusivity<compressible::turbulenceModel>>
        >
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalarField alphap(turbModel.alphaEff(patchI));

    refValue() = massFluxFraction_;
    refGrad() = 0.0;

    valueFraction() =
        1.0
        /
        (
            1.0 +
            alphap*patch().deltaCoeffs()*patch().magSf()/max(mag(phip), SMALL)
        );

    mixedFvPatchField<scalar>::updateCoeffs();

    if (debug)
    {
        scalar phi = gSum(-phip*(*this));

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " mass flux[Kg/s]:" << phi
            << endl;
    }
}

void Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::
setMassFluxFraction(const scalarField& massFluxFrac_org)
{

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        patch().patch()
    );   

    scalarField massFluxFrac = massFluxFrac_org;    
    mpp.distribute(massFluxFrac);


    if (massFluxFrac.size() != massFluxFraction_.size())
    {
        FatalErrorIn
        (
            "Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField"
            "setMassFluxFraction"
            "("
                "const scalarField& "
            ")"
        )
            << " Patch size does not match in the setMassFluxFraction function.\n"
            << nl << nl
            << "There seems to be a serious problem in the case set-up. The patches on the "
            << "default region (region0) and the pyro region should be of the same size." 
            << " Please use a different BC on the patch." << exit(FatalError);
    }
    massFluxFraction_ = massFluxFrac;
}

void Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    massFluxFraction_.writeEntry("massFluxFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField
    );

}

// ************************************************************************* //

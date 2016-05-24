/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "localFuelRadFractionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "zeroGradientFvPatchFields.H"
#include "basicMultiComponentMixture.H"

#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(localFuelRadFractionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            localFuelRadFractionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::localFuelRadFractionEmission::localFuelRadFractionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
//    speciesNames_(0),
//    specieIndex_(label(0)),
//    lookUpTablePtr_(),
//    thermo_(mesh.lookupObject<basicThermo>("thermophysicalProperties")),
//    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
//    Yj_(nSpecies_)
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
    radScaling(coeffsDict_.lookupOrDefault<Switch>("radScaling",false)),
    //fuelStreams_(coeffsDict_.lookup("fuelStreams")),
    //fuelRadFracs_(coeffsDict_.lookup("fuelStreamsRadFracs")),

    fuelStreams_
    (
      IOobject
      (
        "fuelStreams",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      coeffsDict_.lookup("fuelStreams")
    ),

    fuelRadFracs_
    (
      IOobject
      (
        "fuelStreamsRadFracs",
         mesh.time().timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
      ),
      coeffsDict_.lookup("fuelStreamsRadFracs")
    ),

    localRadFrac_
    (
        IOobject
        (
            "localRadFrac",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    )    

{
   
  if (fuelStreams_.size() != fuelRadFracs_.size() )
  {
        FatalErrorIn("localFuelRadFractionEmission absorption/emission model: ")   
                << nl << "Entries in fuelStreams don't match with fuelStreamsRadFracs in "<< coeffsDict_ << " dictionary.." 
                << nl << "Entries in fuelStreams: " << fuelStreams_.size() << nl
                << "Entries in fuelStreamsRadFracs: " << fuelRadFracs_.size()   << nl
                << "Number of entries should be identical " << nl
                << abort(FatalError);

     //Info << "Fatal error " << endl;
  }

  volScalarField ySum(0.*localRadFrac_);
  localRadFrac_ == dimensionedScalar("zero", localRadFrac_.dimensions(), 0.0);

  forAll(fuelStreams_,sID)
  {
    if (mesh.foundObject<volScalarField>(fuelStreams_[sID]))
    {
        const volScalarField& YFuel =
            mesh_.lookupObject<volScalarField>(fuelStreams_[sID]);

        ySum = ySum + YFuel;
        localRadFrac_ = localRadFrac_ + fuelRadFracs_[sID]*YFuel;
    }
  }
  
  localRadFrac_ = localRadFrac_ / max(ySum,1.e-16);

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::localFuelRadFractionEmission::~localFuelRadFractionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::localFuelRadFractionEmission::aCont(const label bandI) const
{

    tmp<volScalarField> a
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            zeroGradientFvPatchVectorField::typeName            
        )
    );

    return a;

}


Foam::tmp<Foam::volScalarField>
Foam::radiation::localFuelRadFractionEmission::eCont(const label bandI) const
{
    tmp<volScalarField> e
    (
        new volScalarField
        (
            IOobject
            (
                "eCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    return e;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::localFuelRadFractionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    //scalar RadFraction = 0;

    //if (radScaling)
    //{
        //const label patch1I = mesh_.boundaryMesh().findPatchID(patchName1_);
        //if(patch1I<0)
        //{
        //    FatalErrorIn("radScaling.H")   
        //        << "patch " << patchName1_ << " not found" << nl
        //        << abort(FatalError);
        //}       
        //const label patch2I = mesh_.boundaryMesh().findPatchID(patchName2_);      
        //if(patch2I<0)
        //{
        //    FatalErrorIn("radScaling.H")   
        //        << "patch " << patchName2_ << " not found" << nl
        //        << abort(FatalError);
        //}       
          
      //  const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");
        //scalar mlr1 = -gSum(phi.boundaryField()[patch1I]);
        //scalar mlr2 = -gSum(phi.boundaryField()[patch2I]);    
        //Info << "mlr for patch " << patchName1_ << " is " << mlr1 << endl;  
        //Info << "mlr for patch " << patchName2_ << " is " << mlr2 << endl;                   
  
        //RadFraction = (mlr1*Ehrr1_ + mlr2*Ehrr2_)
        //            / max(SMALL, (mlr1 + mlr2));        
    //}
    //else
    //{
    //    RadFraction = EhrrCoeff_;
    //}

    if (mesh_.foundObject<volScalarField>("dQ"))
    {
        const volScalarField& dQ =
            mesh_.lookupObject<volScalarField>("dQ");

        if (radScaling)
        {
          E.ref().ref() = localRadFrac_*dQ;
        }
        else
        {
          E.ref().ref() = EhrrCoeff_*dQ;
        }

    }

    return E;
}


// ************************************************************************* //

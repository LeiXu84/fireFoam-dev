/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

#include "eddyDissMultiFuelStreamModel.H"
#include "scalarIOList.H"
#include "wordIOList.H"

namespace Foam
{
namespace combustionModels
{

// * * * Initialization routine * * *

template<class CombThermoType, class ThermoType>
void eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::initPars()
{
 
  if (this->semiImplicit_)
  {
        FatalErrorIn
        (
            "eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::"
            "eddyDissMultiFuelStreamModel"
            "("
                "const word&, "
                "const fvMesh&"
            ")"
        )
            << " semiImplicit_ option currently not supported for eddyDissMultiFuelStreamModel model:\n"
            << nl << nl
            << "Please set semiImplicit_ option to no in the " << this->coeffs().name() 
            << " dictionary entry" << exit(FatalError);

  }

  forAll(fuelSpecies_,sID)
  {
      fuelSpeciesID_[sID] = this->singleMixturePtr_->species()[fuelSpecies_[sID]];
  } 
  //Info << "fuel species names: " << fuelSpecies_ << endl;
  //Info << "fuel species ID: " << fuelSpeciesID_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::eddyDissMultiFuelStreamModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& combustionProperties,
    const word& phaseName
)
:
    singleStepCombustion<CombThermoType, ThermoType>
    (
        modelType,
        mesh,
        combustionProperties,
        phaseName
    ),
    C_(readScalar(this->coeffs().lookup("C_EDC"))),
    Cd_(readScalar(this->coeffs().lookup("C_Diff"))),
    Cstiff_(readScalar(this->coeffs().lookup("C_Stiff"))),
    fuelSpecies_(this->coeffs().lookup("fuelSpecies")),
    fuelSpeciesID_(fuelSpecies_.size()),

    fuelSpeciesSum_
    (
        IOobject
        (
            "fuelSpeciesSum",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    )
{
    initPars();
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::~eddyDissMultiFuelStreamModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::rtTurb() const
{
    return C_*this->turbulence().epsilon()/
              max(this->turbulence().k(),
              dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));
}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::rtDiff() const
{
    const volScalarField& YO2 = this->thermoPtr_->composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
         turbulenceModel::propertiesName
        );

    return Cd_*this->thermoPtr_->alpha()/this->rho()/sqr(lesModel.delta());
}

template<class CombThermoType, class ThermoType>
void eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {


        this->singleMixturePtr_->fresCorrect();

        //const label fuelI = this->singleMixturePtr_->fuelIndex();

        //const volScalarField& YFuel =
        //    this->thermoPtr_->composition().Y()[fuelI];

        fuelSpeciesSum_ == dimensionedScalar("zero", dimless, 0.0);
        volScalarField& YFuel = fuelSpeciesSum_; // 0.*this->thermoPtr_->composition().Y()[fuelI];
        //volScalarField YFuel = 0.*this->thermoPtr_->composition().Y()[fuelI];

        forAll(fuelSpecies_,sID)
        {
          YFuel == YFuel + this->thermoPtr_->composition().Y()[fuelSpeciesID_[sID]];
        }

        //const volScalarField& YFuel =
        //    this->thermoPtr_->composition().Y()[fuelSpeciesID_[0]];

        //fuelSpeciesSum_ = YFuel;

        // Code to set the appropriate localRadFrac to be used in the localFuelRadFractionEmission model..
       if (fuelSpeciesSum_.mesh().foundObject<wordIOList>("fuelStreams"))
       //if (this->mesh().template foundObject<wordList>("fuelStreams"))
       {
          const  wordIOList& fStreams =
             fuelSpeciesSum_.mesh().lookupObject<wordIOList>("fuelStreams");
        //     this->mesh().template lookupObject<wordList>("fuelStreams");

          const scalarIOList& fStreamsRadFracs =
                fuelSpeciesSum_.mesh().lookupObject<scalarIOList>("fuelStreamsRadFracs");
          //  this->mesh().template lookupObject<scalarList>("fuelStreamsRadFracs");

          volScalarField& localRadFrac =
             const_cast<volScalarField&>(fuelSpeciesSum_.mesh().lookupObject<volScalarField>("localRadFrac"));
          //  const_cast<volScalarField&>(this->mesh().template lookupObject<volScalarField>("localRadFrac"));

          localRadFrac == dimensionedScalar("zero", localRadFrac.dimensions(), 0.0);


          forAll(fStreams,fsID)
          { 
            forAll(fuelSpecies_,sID)
            {
              if (fuelSpecies_[sID] == fStreams[fsID])
              {
                localRadFrac = localRadFrac + fStreamsRadFracs[fsID]*this->thermoPtr_->composition().Y()[fuelSpeciesID_[sID]]; 
              }
            }  
            
          }

          localRadFrac = localRadFrac / max(YFuel,1.e-299);
          localRadFrac = max(0.,localRadFrac);
       }

        const dimensionedScalar s = this->singleMixturePtr_->s();

        if (this->thermoPtr_->composition().contains("O2"))
        {
            const volScalarField& YO2 = this->thermoPtr_->composition().Y("O2");

//            this->wFuel_ ==
//                this->rho()/(this->mesh().time().deltaT()*C_)
//               *min(YFuel, YO2/s.value());

/*
            this->wFuel_ ==
                  C_
                * this->rho()
                * this->turbulence().epsilon()
                / max(this->turbulence().k(),
                  dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL))
                * min(YFuel, YO2/s.value());
*/

/*
            this->wFuel_ ==
                  this->rho()
                * min(YFuel, YO2/s.value())
                * max(rtTurb(),rtDiff());
*/

            volScalarField rt = max(rtTurb(),rtDiff()); 

/*
            this->wFuel_ ==
                  this->rho()
                * min(YFuel, YO2/s.value())
                / this->mesh_.time().deltaT() / Cstiff_
                * (1 - exp(- Cstiff_*this->mesh_.time().deltaT() * rt));
*/

/*
            this->wFuel_ ==
                  min(this->rho()
                * min(max(0*YFuel,YFuel), max(0*YO2,YO2)/s.value())
                / this->mesh_.time().deltaT() / Cstiff_
                * (1 - exp(- Cstiff_*this->mesh_.time().deltaT() * rt)),
                this->rho()*min(max(0*YFuel,YFuel), max(0*YO2,YO2)/s.value())/this->mesh().time().deltaT());
*/

// LINES BELOW ARE THE ACTUAL CODE.. ONLY COMMENTING HERE FOR DEBUGGING PURPOSES
            this->wFuel_ ==
                   this->rho()
                * min(max(0*YFuel,YFuel), max(0*YO2,YO2)/s.value())
                / this->mesh_.time().deltaT() *
                min(1./ Cstiff_* (1 - exp(- Cstiff_*this->mesh_.time().deltaT() * rt)),1.0);

        }
    }
}

//

template<class CombThermoType, class ThermoType>
tmp<fvScalarMatrix> eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::R
(
    volScalarField& Y
) const
{
    const label specieI = this->thermoPtr_->composition().species()[Y.name()];

    volScalarField wSpecie
    (
        this->wFuel_*this->singleMixturePtr_->specieStoichCoeffs()[specieI]
    );

    const label fuelI = this->singleMixturePtr_->fuelIndex();

    if (specieI == fuelI)
    {
      wSpecie == 0.*this->wFuel_;
    }

    forAll(fuelSpeciesID_,specID)
    {
      if (fuelSpeciesID_[specID] == specieI)
      {

/*
        forAll(fuelSpeciesSum_.internalField(),cellI)
        {
          if (fuelSpeciesSum_.internalField()[cellI] <= 0.) 
          {
              //wSpecie.internalField()[cellI] = 0.;
          }
          else
          {
              scalar testRat = Y.internalField()[cellI] / fuelSpeciesSum_.internalField()[cellI];
              wSpecie.internalField()[cellI] = this->wFuel_.internalField()[cellI] * this->singleMixturePtr_->specieStoichCoeffs()[fuelI] *
                                                Y.internalField()[cellI] / fuelSpeciesSum_.internalField()[cellI]; 
                                               //testRat; // (1.+1.e-16); // Y.internalField()[cellI] / fuelSpeciesSum_.internalField()[cellI]; 
              Info << cellI << " " << Y.internalField()[cellI] / fuelSpeciesSum_.internalField()[cellI] << endl;
          } 
        }
*/

        volScalarField testRat = Y / max(1.e-299,fuelSpeciesSum_);        
        //wSpecie == this->wFuel_ * this->singleMixturePtr_->specieStoichCoeffs()[fuelI] * Y / max(1.e-299,fuelSpeciesSum_);
        wSpecie == this->wFuel_ * this->singleMixturePtr_->specieStoichCoeffs()[fuelI] * testRat; 
        // Line below ONLY FOR DEBUGGING.. needs to be corrected later..
        //wSpecie == this->wFuel_ * this->singleMixturePtr_->specieStoichCoeffs()[fuelI]; //  * Y / max(1e-16,fuelSpeciesSum_);

/*
        Info << "InternalField" << endl; 
        forAll(fuelSpeciesSum_.internalField(),cellI)
        {
           Info << wSpecie.internalField()[cellI] << "  " << this->singleMixturePtr_->specieStoichCoeffs()[fuelI]*this->wFuel_.internalField()[cellI] 
           << "  " << wSpecie.internalField()[cellI] - this->singleMixturePtr_->specieStoichCoeffs()[fuelI]*this->wFuel_.internalField()[cellI] << endl;
        }  

        Info << wSpecie - this->wFuel_ * this->singleMixturePtr_->specieStoichCoeffs()[fuelI] << endl;
*/
      }
     
       //fuelSpeciesSum_ == wSpecie - this->wFuel_ * this->singleMixturePtr_->specieStoichCoeffs()[fuelI];
       //Info << wSpecie << endl;
       //Info << this->wFuel_ * this->singleMixturePtr_->specieStoichCoeffs()[fuelI] << endl;    
 
    }

    if (this->semiImplicit_)
    {
        const label fNorm = this->singleMixturePtr_->specieProd()[specieI];
        const volScalarField fres(this->singleMixturePtr_->fres(specieI));
        wSpecie /= max(fNorm*(Y - fres), scalar(1e-2));

        return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
    }
    else
    {
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
}

// 

template<class CombThermoType, class ThermoType>
tmp<volScalarField>
eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::Sh() const
{
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermoPtr_->composition().Y(fuelI));
    tmp<volScalarField> shSum(0.*this->singleMixturePtr_->qFuel()*(R(YFuel) & YFuel));

    forAll(fuelSpeciesID_,specID)
    {
        volScalarField& YFuel1 =
          const_cast<volScalarField&>(this->thermoPtr_->composition().Y(fuelSpeciesID_[specID]));

        shSum = shSum -this->singleMixturePtr_->qFuel()*(R(YFuel1) & YFuel1);
        Info << " qFuel: " << this->singleMixturePtr_->qFuel() << endl;
     }
  
    return shSum;
    //return this->wFuel_*this->singleMixturePtr_->qFuel();
}

//

template<class CombThermoType, class ThermoType>
tmp<volScalarField>
eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
        )
    );

    if (this->active())
    {
        volScalarField& dQ = tdQ.ref();
        dQ.ref() = this->mesh().V()*Sh()();
    }
    return tdQ;
}

//

template<class CombThermoType, class ThermoType>
bool eddyDissMultiFuelStreamModel<CombThermoType, ThermoType>::read()
{
    if (singleStepCombustion<CombThermoType, ThermoType>::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

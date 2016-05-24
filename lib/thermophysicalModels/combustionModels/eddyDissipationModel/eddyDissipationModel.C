/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "eddyDissipationModel.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationModel<CombThermoType, ThermoType>::eddyDissipationModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& phaseName
)
:
    singleStepCombustion<CombThermoType, ThermoType>
    (
        modelType,
        mesh,
        phaseName
    ),
    C_(readScalar(this->coeffs().lookup("C_EDC"))),
    Cd_(readScalar(this->coeffs().lookup("C_Diff"))),
    Cstiff_(readScalar(this->coeffs().lookup("C_Stiff")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationModel<CombThermoType, ThermoType>::~eddyDissipationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationModel<CombThermoType, ThermoType>::rtTurb() const
{
    return C_*this->turbulence().epsilon()/
              max(this->turbulence().k(),
              dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));
}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationModel<CombThermoType, ThermoType>::rtDiff() const
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
void eddyDissipationModel<CombThermoType, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {
        this->singleMixturePtr_->fresCorrect();

        const label fuelI = this->singleMixturePtr_->fuelIndex();

        const volScalarField& YFuel = this->thermoPtr_->composition().Y()[fuelI];

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

            volScalarField rt(max(rtTurb(),rtDiff())); 

            // clipping of wFuel to prevent negative HRR
            // this->wFuel_ ==
            //     this->rho()
            //     * min(max(0*YFuel,YFuel), max(0*YO2,YO2)/s.value())
            //     / this->mesh_.time().deltaT() *
            //     min(1./ Cstiff_* (1 - exp(- Cstiff_*this->mesh_.time().deltaT() * rt)),1.0);
            
            this->wFuel_ ==
                  this->rho()
                * min(YFuel, YO2/s.value())
                / this->mesh_.time().deltaT() / Cstiff_
                * (1 - exp(- Cstiff_*this->mesh_.time().deltaT() * rt));
        }
    }
}


template<class CombThermoType, class ThermoType>
bool eddyDissipationModel<CombThermoType, ThermoType>::read()
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

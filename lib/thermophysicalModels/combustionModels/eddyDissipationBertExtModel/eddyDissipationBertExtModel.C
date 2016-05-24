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

#include "eddyDissipationBertExtModel.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "fvCFD.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationBertExtModel<CombThermoType, ThermoType>::eddyDissipationBertExtModel
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
    Cstiff_(readScalar(this->coeffs().lookup("C_Stiff"))),
    tExt_(this->coeffs().lookupOrDefault("ExtinctionStart", 5.0)),
    TFuelExt_(this->coeffs().lookupOrDefault("FuelExtTemp", 400.0)),
    TFuelStarExt_(this->coeffs().lookupOrDefault("FuelStarExtTemp", 1000.0)),
    Cstrain_(this->coeffs().lookupOrDefault("Cstrain", 0.25)),
    Cevap_(this->coeffs().lookupOrDefault("Cevap", 0.5)),
    dQMin_(this->coeffs().lookupOrDefault("dQMin", 1.0e6)),
    ignitionLocation_(this->coeffs().lookupOrDefault("ignitionLocation", 0.5)),
    flameQ_
    (
        IOobject
        (
            "flameQ",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
	zeroGradientFvPatchScalarField::typeName
    ),
    FEF_
    (
        IOobject
        (
            "FEF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    FIF_
    (
        IOobject
        (
            "FIF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 1.0)
    ),
    strainRate_
    (
        IOobject
        (
            "strainRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    ),
    Textinction_
    (
        IOobject
        (
            "Textinction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature, 0.0)
    ),
    hsAdditional_
    (
        IOobject
        (
            "hsAdditional",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimMass/dimVolume, 0.0)
    ),
    entrainRate_
    (
        IOobject
        (
            "entrainRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    ),
    WCO2_
    (
        IOobject
        (
            "WCO2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    WO2_
    (
        IOobject
        (
            "WO2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    WH2O_
    (
        IOobject
        (
            "WH2O",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    WF_
    (
        IOobject
        (
            "WF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    WFstar_
    (
        IOobject
        (
            "WFstar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    WFNet_
    (
        IOobject
        (
            "WFNet",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    WFstarNet_
    (
        IOobject
        (
            "WFstarNet",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    evaporationRate_
    (
        IOobject
        (
            "evaporationRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationBertExtModel<CombThermoType, ThermoType>::~eddyDissipationBertExtModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationBertExtModel<CombThermoType, ThermoType>::rtTurb() const
{
    return C_*this->turbulence().epsilon()/
              max(this->turbulence().k(),
              dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));
}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationBertExtModel<CombThermoType, ThermoType>::rtDiff() const
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
void eddyDissipationBertExtModel<CombThermoType, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {
        this->singleMixturePtr_->fresCorrect();

        const label fuelI = this->singleMixturePtr_->fuelIndex();

        const volScalarField& YFuel =
            this->thermoPtr_->composition().Y()[fuelI];

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

            this->wFuel_ ==
                  this->rho()
                * min(YFuel, YO2/s.value())
                / this->mesh_.time().deltaT() / Cstiff_
                * (1 - exp(- Cstiff_*this->mesh_.time().deltaT() * rt));

	    //- Flame Extinction using Bert's Model
    	    //
    	    //- Strain Rate
    	    strainRate_ = Cstrain_*rt;
    	    
    	    //- Species Index
    	    label indexH2O(this->thermoPtr_->composition().species()["H2O"]);
    	    label indexCO2(this->thermoPtr_->composition().species()["CO2"]);
    	    label indexN2(this->thermoPtr_->composition().species()["N2"]);
    	    label indexO2(this->thermoPtr_->composition().species()["O2"]);
    	    const volScalarField& YH2O = this->thermoPtr_->composition().Y("H2O");
    	    const volScalarField& YCO2 = this->thermoPtr_->composition().Y("CO2");
    	    const volScalarField& YN2 = this->thermoPtr_->composition().Y("N2");
    	    const volScalarField& YFstar = this->thermoPtr_->composition().Y("Fstar");
    	    	
    	    //- Reaction rates
    	    WF_ = this->wFuel_;
    	    WCO2_ = this->singleMixturePtr_->specieStoichCoeffs()[indexCO2]*this->wFuel_;
    	    WH2O_ = this->singleMixturePtr_->specieStoichCoeffs()[indexH2O]*this->wFuel_;
    	    WO2_ = this->singleMixturePtr_->specieStoichCoeffs()[indexO2]*this->wFuel_;
    	    dimensionedScalar qF(this->singleMixturePtr_->qFuel());
    	    //scalar qF(this->singleMixturePtr_->qFuel());
    	    
    	    //- Spray Info
    	    const volScalarField sprayDensity = 
    	    	this->mesh().template lookupObject<volScalarField>("rhoSpray");
    	    const volScalarField dQFuel = 
    	    	this->mesh().template lookupObject<volScalarField>("dQ");
    	    const volScalarField TCell = this->thermoPtr_->T();
    	    const volScalarField pCell = this->thermoPtr_->p();
    	    const volScalarField rhoCell = this->thermoPtr_->rho();
    	    //const volScalarField& rhoCellRef = this->thermoPtr_->rho();
    	    //const volScalarField& TCellRef = this->thermoPtr_->T();
    	    //const volScalarField& pCellRef = this->thermoPtr_->p();
            scalar Hwater(2.257e6);

	    //- Transport equation for flame surface
	    const compressible::LESModel& lesModel =
            YO2.db().lookupObject<compressible::LESModel>
	    (
	        turbulenceModel::propertiesName
	    );
	    //const surfaceScalarField alphaRhoPhi = lesModel.alphaRhoPhi();
	    //const surfaceScalarField PhiLES = lesModel.phi();
	    const volScalarField muEffective = lesModel.muEff();
    	    const surfaceScalarField phi = 
    	    	this->mesh().template lookupObject<surfaceScalarField>("phi");
	    flameQ_ = dQFuel;
	    flameQ_ = flameQ_ + (
		    fvc::laplacian(muEffective, flameQ_) - fvc::div(phi, flameQ_)
		    )*this->mesh_.time().deltaT()/rhoCell;
	    const vectorField& cellCenter = this->mesh_.cellCentres();
    	    
    	    //Textinction_ = 1368.52*pow(strainRate_, 0.1055);
    	    forAll(Textinction_,cellI)
    	    {
		//- For CH4
    	    	//if(strainRate_[cellI] > 10)
    	    	//{
    	    	//    Textinction_[cellI] = 1368.52*pow(strainRate_[cellI], 0.1055);
    	    	//}
    	    	//else
    	    	//{
    	    	//    Textinction_[cellI] = 1368.52*pow(10, 0.1055);
    	    	//}
		
		//- For C3H8
    	    	if(strainRate_[cellI] > 7)
    	    	{
    	    	    Textinction_[cellI] = 1328.5*pow(strainRate_[cellI], 0.1143);
    	    	}
    	    	else
    	    	{
    	    	    //Textinction_[cellI] = 1328.5*pow(7, 0.1143);
    	    	    Textinction_[cellI] = 1659;
    	    	}
    	    
    	    	entrainRate_[cellI] = -WO2_[cellI]/(YO2[cellI]*rhoCell[cellI]+SMALL);

		scalar pC(pCell[cellI]);
		scalar TC(TCell[cellI]);
		scalar TEC(Textinction_[cellI]);
    	    
    	    	scalar hH2OText = 
    	    	    this->thermoPtr_->composition().Hs(indexH2O,pC,TEC);
    	    	scalar hH2OTcur = 
    	    	    this->thermoPtr_->composition().Hs(indexH2O,pC,TC);
    	    	scalar hCO2Text = 
    	    	    this->thermoPtr_->composition().Hs(indexCO2,pC,TEC);
    	    	scalar hCO2Tcur = 
    	    	    this->thermoPtr_->composition().Hs(indexCO2,pC,TC);
    	    	scalar hN2Text = 
    	    	    this->thermoPtr_->composition().Hs(indexN2,pC,TEC);
    	    	scalar hN2Tcur = 
    	    	    this->thermoPtr_->composition().Hs(indexN2,pC,TC);

    	    	evaporationRate_[cellI] = Cevap_*entrainRate_[cellI] * sprayDensity[cellI];
    	    
    	    	hsAdditional_[cellI] = evaporationRate_[cellI]*(Hwater + hH2OText - hH2OTcur)
    	    	    + entrainRate_[cellI]*rhoCell[cellI]*
    	    	    (YN2[cellI]*(hN2Text-hN2Tcur)
    	    	    + YH2O[cellI]*(hH2OText-hH2OTcur)
    	    	    + YCO2[cellI]*(hCO2Text-hCO2Tcur))
    	    	    + WCO2_[cellI]*(hCO2Text-hCO2Tcur) + WH2O_[cellI]*(hH2OText-hH2OTcur)
    	    	    //- dQFuel[cellI];
    	    	    //- this->Sh()()[cellI];//- A bug here, this is the actual HRR, not from wFuel_;
    	    	    - WF_[cellI]*qF.value();
    	    
    	    	    ////- Assume Reignition is always harder than extinction
    	    	    if(TCell[cellI] > TFuelStarExt_)
    	    	    {
    	    	        FIF_[cellI] = 1;
    	    	        FEF_[cellI] = 0;
    	    	    }
    	    	    else	//- might be a bug: if no extinction, re-ignition occurs?
    	    	    {
    	    	        FIF_[cellI] = 0;
    	    	        //if(hsAdditional_[cellI] > 1)
		        if((this->mesh_.time().value() > tExt_)
		            && (TCell[cellI] < TFuelExt_ || hsAdditional_[cellI] > 1))
    	    	        {
    	    	        	FEF_[cellI] = 1;
    	    	        }
    	    	        else
    	    	        {
    	    	        	FEF_[cellI] = 0;
    	    	        }
    	    	    }

		    //- New approach based on flame surface transport
		    //if(flameQ_[cellI] > dQMin_)
		    //{
		    //    FIF_[cellI] = 1;
		    //    FEF_[cellI] = 0;
		    //}
		    //else
		    //{
		    //    FIF_[cellI] = 0;
		    //    FEF_[cellI] = 1;
		    //}
		    //if(cellCenter[cellI].z() > ignitionLocation_)
		    //{
		    //    FIF_[cellI] = 1;
		    //    FEF_[cellI] = 0;
		    //}
    	    }
                this->WFstar_ ==
                  this->rho()*YFstar
                / this->mesh_.time().deltaT() / Cstiff_
                * (1 - exp(- Cstiff_*this->mesh_.time().deltaT() * rt));

	    //this->WFNet_ == R(YFuel) & YFuel;
	    //this->WFstarNet_ == R(YFstar) & YFstar;
	    this->WFNet_ == -this->wFuel_ + this->WFstar_*FIF_;
	    this->WFstarNet_ == this->wFuel_*FEF_ - this->WFstar_*FIF_;
        }
    }
}


template<class CombThermoType, class ThermoType>                                                    
tmp<fvScalarMatrix>  
eddyDissipationBertExtModel<CombThermoType, ThermoType>::R
(
    volScalarField& Y
) const
{
    const label specieI = this->thermoPtr_->composition().species()[Y.name()]; 
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    if(specieI == fuelI)
    {
    	volScalarField wSpecie
    	(
	    this->wFuel_*this->singleMixturePtr_->specieStoichCoeffs()[specieI]
		+ FIF_*this->WFstar_
    	);
	return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
    else if(Y.name() == "Fstar")
    {
    	volScalarField wSpecie
    	(
	   FEF_*this->wFuel_ - FIF_*this->WFstar_
    	);
	return wSpecie + fvm::Sp(0.0*wSpecie, Y); 
    }
    else
    {
    	volScalarField wSpecie
    	(
	   (1-FEF_)*this->wFuel_*this->singleMixturePtr_->specieStoichCoeffs()[specieI]
    	);
	return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
}                               


template<class CombThermoType, class ThermoType>
tmp<volScalarField>                     
eddyDissipationBertExtModel<CombThermoType, ThermoType>::Sh() const
{                                       
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    volScalarField& YFuel =             
        const_cast<volScalarField&>(this->thermoPtr_->composition().Y(fuelI));
	label indexFstar(this->thermoPtr_->composition().species()["Fstar"]);
    volScalarField& YFstar =             
        const_cast<volScalarField&>(this->thermoPtr_->composition().Y(indexFstar));
    return -this->singleMixturePtr_->qFuel()*((R(YFuel) & YFuel) + (R(YFstar) & YFstar));
}


template<class CombThermoType, class ThermoType>
bool eddyDissipationBertExtModel<CombThermoType, ThermoType>::read()
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

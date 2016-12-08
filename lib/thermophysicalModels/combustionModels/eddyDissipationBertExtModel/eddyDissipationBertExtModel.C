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
    tExt_(this->coeffs().lookupOrDefault("ExtinctionStart", 5.0)),
    TFuelExt_(this->coeffs().lookupOrDefault("FuelExtTemp", 400.0)),
    TFuelStarExt_(this->coeffs().lookupOrDefault("FuelStarExtTemp", 1000.0)),
    Cstrain_(this->coeffs().lookupOrDefault("Cstrain", 0.25)),
    Cevap_(this->coeffs().lookupOrDefault("Cevap", 0.5)),
    XrExtinction_(this->coeffs().lookupOrDefault("XrExtinction", 0.05)),
    nearWallExtinction_(this->coeffs().lookupOrDefault("nearWallExtinction", false)),
    radiativeHeatLoss_(this->coeffs().lookupOrDefault("radiativeHeatLoss", false)),
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
    Tad_
    (
        IOobject
        (
            "Tad",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature, 293.0)
    ),
    TadLoss_
    (
        IOobject
        (
            "TadLoss",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature, 293.0)
    ),
    WCO2_
    (
        IOobject
        (
            "WCO2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    nearWallCells_
    (
        IOobject
        (
            "nearWallCells",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    )
{
    const polyBoundaryMesh& bm = mesh.boundaryMesh();
    const wordList patchNames = bm.names();
    const wordList patchTypes = bm.types();
    forAll(patchNames,nameI)
    {
    label patchIndex = bm.findPatchID(patchNames[nameI]);
    Info<<"Name I = "<<patchIndex<<tab<<patchNames[nameI]<<tab<<patchTypes[nameI]<<endl;
    if((patchTypes[nameI] == "mappedWall") || (patchNames[nameI] == "burner"))
    {
        Info<<"Near wall cells for: "<<patchNames[nameI]<<tab<<patchTypes[nameI]<<endl;
        const polyPatch& pp = bm[patchIndex];
            forAll(pp, faceI)
            {
                label globalFaceI = pp.start() + faceI;
                label ownCellI = mesh.faceOwner()[globalFaceI];
                nearWallCells_[ownCellI] = 1.0;
            }
    }
    }

}


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

            calculateFlameTemperature();
            
            //- Species Index
            label indexH2O(this->thermoPtr_->composition().species()["H2O"]);
            label indexCO2(this->thermoPtr_->composition().species()["CO2"]);
            label indexN2(this->thermoPtr_->composition().species()["N2"]);
            label indexO2(this->thermoPtr_->composition().species()["O2"]);
            const volScalarField& YH2O = this->thermoPtr_->composition().Y("H2O");
            const volScalarField& YCO2 = this->thermoPtr_->composition().Y("CO2");
            const volScalarField& YN2 = this->thermoPtr_->composition().Y("N2");
            const volScalarField& YFstar = this->thermoPtr_->composition().Y("Fstar");

            const volScalarField& TCellRef = this->thermoPtr_->T();
                
            //- Reaction rates
            WF_ = this->wFuel_;
            WCO2_ = this->singleMixturePtr_->specieStoichCoeffs()[indexCO2]*this->wFuel_;
            WH2O_ = this->singleMixturePtr_->specieStoichCoeffs()[indexH2O]*this->wFuel_;
            WO2_ = this->singleMixturePtr_->specieStoichCoeffs()[indexO2]*this->wFuel_;
            dimensionedScalar qF(this->singleMixturePtr_->qFuel());

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
                    Textinction_[cellI] = 1659;
                }

        //- Re-ignition
        if(TCellRef[cellI] > TFuelStarExt_)
        {
            FIF_[cellI] = 1;
        }
        else
        {
            FIF_[cellI] = 0;
        }
            
        //- Extinction
        scalar TcellI = Tad_[cellI];
        if(radiativeHeatLoss_)
        {
            TcellI = TadLoss_[cellI];
        }
        if( 
            (this->mesh_.time().value() > tExt_)
            && 
            (TCellRef[cellI] < TFuelExt_ || TcellI < Textinction_[cellI])
          )
        {
            FEF_[cellI] = 1;

            if(!nearWallExtinction_)
            {
            if(nearWallCells_[cellI] == 1.0)
            {
                FEF_[cellI] = 0;
            }
            }

            if((YFuel[cellI] > 1e-3) && (YO2[cellI] > 1e-3))
            {
            FIF_[cellI] = 0;
            }
        }
        else
        {
            FEF_[cellI] = 0;
            //FIF_[cellI] = 1;
        }
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

    const label indexFstar(this->thermoPtr_->composition().species()["Fstar"]);
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

template<class CombThermoType, class ThermoType>
void eddyDissipationBertExtModel<CombThermoType, ThermoType>::calculateFlameTemperature()
{
    Info<<"Calculate Adaiabatic Flame Temperature!"<<endl;
    //- Get species index
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    label indexH2O(this->thermoPtr_->composition().species()["H2O"]);
    label indexCO2(this->thermoPtr_->composition().species()["CO2"]);
    label indexN2(this->thermoPtr_->composition().species()["N2"]);
    label indexO2(this->thermoPtr_->composition().species()["O2"]);
    label indexFstar(this->thermoPtr_->composition().species()["Fstar"]);

    //- Get species mass fraction
    const volScalarField& YFuel = this->thermoPtr_->composition().Y()[fuelI];
    const volScalarField& YO2 = this->thermoPtr_->composition().Y("O2");
    const volScalarField& YN2 = this->thermoPtr_->composition().Y("N2");
    const volScalarField& YCO2 = this->thermoPtr_->composition().Y("CO2");
    const volScalarField& YH2O = this->thermoPtr_->composition().Y("H2O");
    const volScalarField& YFstar = this->thermoPtr_->composition().Y("Fstar");

    const volScalarField& TCellRef = this->thermoPtr_->T();
    const volScalarField& pCellRef = this->thermoPtr_->p();
    dimensionedScalar qF(this->singleMixturePtr_->qFuel());
    const dimensionedScalar s = this->singleMixturePtr_->s();
    const volScalarField rhoCell = this->thermoPtr_->rho();

    //- Spray Info
    const volScalarField sprayDensity = 
        this->mesh().template lookupObject<volScalarField>("rhoSpray");

    //- Hard coded, for C3H8 only.
    forAll(YFuel, cellI)
    {
    if((YFuel[cellI] > 1e-3) && (YO2[cellI] > 1e-3))
    {
        //scalar mCell(this->rho()[cellI]);
        scalar mCell(rhoCell[cellI]);
        scalar mFuel(mCell*(YFuel[cellI] + YFstar[cellI]));
        scalar mO2(mCell*YO2[cellI]);
        scalar mN2(mCell*YN2[cellI]);
        scalar mCO2(mCell*YCO2[cellI]);
        scalar mH2O(mCell*YH2O[cellI]);
        scalar mWater(sprayDensity[cellI]);

        scalar nFuel(mFuel/44.0);
        scalar nO2(mO2/32.0);
        scalar nN2(mN2/28.0);
        scalar nCO2(mCO2/44.0);
        scalar nH2O(mH2O/18.0);

        scalar RN(nN2*5.0/(nO2+1.6667*nCO2));
        scalar nN2r(RN*nO2/5.0);
        scalar nN2p(max(0.0, (nN2 - nN2r)));
        scalar nH2Op(1.33333*nCO2);
        scalar nH2Ow(max(0.0, (nH2O - nH2Op)));
        scalar mN2r(nN2r*28.0);
        scalar mN2p(nN2p*28.0);
        scalar mH2Op(nH2Op*18.0);
        scalar mH2Ow(nH2Ow*18.0);

        scalar dmFuel(min(mFuel, mO2/s.value()));
        scalar dmO2(dmFuel*s.value());
        scalar dmCO2(3.0*dmFuel);
        scalar dmH2Op(4.0*dmFuel*18.0/44.0);
        scalar dmN2(RN*dmFuel*28.0/44.0);

        scalar dnFuel(dmFuel/44.0);
        scalar entrainRatio(dnFuel*(1+5.0+RN)/(nO2+nN2r+nFuel));
        scalar dmWater(Cevap_*mWater*entrainRatio);

        //- Coefficients for the equation
        scalar CCO2(dmCO2 + entrainRatio*mCO2);
        scalar CH2O(dmH2Op + entrainRatio*mH2O + dmWater);
        scalar CN2(dmN2 + entrainRatio*mN2p);
        scalar hrr(dmFuel*qF.value());
        scalar hrrLoss(dmFuel*qF.value()*(1.0-XrExtinction_));

        //aN2_[cellI] = RN;
        //entrainRatio_[cellI] = entrainRatio;
        //dmCO2_[cellI] = dmCO2;
        //mCO2_[cellI] = mCO2;
        //dmH2O_[cellI] = dmH2Op;
        //mH2O_[cellI] = mH2O;
        //dmN2_[cellI] = dmN2;
        //mN2_[cellI] = mN2p;

            scalar hF = 
                this->thermoPtr_->composition().Hs(fuelI,pCellRef[cellI],TCellRef[cellI])
        -
                this->thermoPtr_->composition().Hs(fuelI,pCellRef[cellI],293.15);
            scalar hO2 = 
                this->thermoPtr_->composition().Hs(indexO2,pCellRef[cellI],TCellRef[cellI])
        -
                this->thermoPtr_->composition().Hs(indexO2,pCellRef[cellI],293.15);
            scalar hCO2 = 
                this->thermoPtr_->composition().Hs(indexCO2,pCellRef[cellI],TCellRef[cellI])
        -
                this->thermoPtr_->composition().Hs(indexCO2,pCellRef[cellI],293.15);
            scalar hH2O = 
                this->thermoPtr_->composition().Hs(indexH2O,pCellRef[cellI],TCellRef[cellI])
        -
                this->thermoPtr_->composition().Hs(indexH2O,pCellRef[cellI],293.15);
            scalar hN2 = 
                this->thermoPtr_->composition().Hs(indexN2,pCellRef[cellI],TCellRef[cellI])
        -
                this->thermoPtr_->composition().Hs(indexN2,pCellRef[cellI],293.15);

        scalar mhp(entrainRatio*(mFuel*hF+mO2*hO2+mCO2*hCO2+mH2O*hH2O+mN2*hN2));
        scalar hEvap(dmWater*2.6e6);

        //mhp_[cellI] = mhp;

        //- Solve flame temperature directly
        //- eA*X^2 + eB*X - eC = eD
        scalar eA(CCO2*0.0926897 + CH2O*0.282627 + CN2*0.0672494);
        scalar eB(CCO2*1077.18 + CH2O*1892.55 + CN2*1064.99);
        scalar eC(CCO2*38985.8 + CH2O*24301.8 + CN2*14618.6);
        scalar eD(mhp + hrr - hEvap);
        scalar eDLoss(mhp + hrrLoss - hEvap);
        //Tad_[cellI] = (sqrt(eB*eB + 4.0*eA*(eC+eD)) - eB)/(2.0*eA) + 293.15;
        TadLoss_[cellI] = (sqrt(eB*eB + 4.0*eA*(eC+eDLoss)) - eB)/(2.0*eA) + 293.15;

        //- Solve flame temperature by iteration
        scalar T1(1500);
        scalar T2(1600);
        scalar T3(1700);
        scalar dH1(0);
        scalar dH2(0);
        scalar dH3(1e6);
        scalar nIter(0);
        dH2 = CCO2*this->thermoPtr_->composition().Hs(indexCO2,pCellRef[cellI],T2)
             + CH2O*this->thermoPtr_->composition().Hs(indexH2O,pCellRef[cellI],T2)
         + CN2*this->thermoPtr_->composition().Hs(indexN2,pCellRef[cellI],T2)
         - eD;
        while((dH3 > CCO2*1e4) && (T3 > 300) && (T3 < 3000) && (nIter < 5))
        {
            dH3 = CCO2*this->thermoPtr_->composition().Hs(indexCO2,pCellRef[cellI],T3)
                 + CH2O*this->thermoPtr_->composition().Hs(indexH2O,pCellRef[cellI],T3)
                 + CN2*this->thermoPtr_->composition().Hs(indexN2,pCellRef[cellI],T3)
                 - eD;
            T1 = T2;
            dH1 = dH2;
            T2 = T3;
            dH2 = dH3;
            T3 = T1 - (T2-T1)*dH1/(dH2-dH1);
            nIter = nIter + 1;
        }
        Tad_[cellI] = T3;
    }
    else
    {
        Tad_[cellI] = 0;
        TadLoss_[cellI] = 0;
    }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

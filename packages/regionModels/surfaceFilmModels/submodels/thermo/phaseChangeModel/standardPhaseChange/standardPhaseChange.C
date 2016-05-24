/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "standardPhaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoSingleLayer.H"
#include "specie.H"
#include "heatTransferModel.H"
#include "filmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardPhaseChange, 0);

addToRunTimeSelectionTable
(
    phaseChangeModel,
    standardPhaseChange,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar standardPhaseChange::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    if (Sc < 0.01)
    {
        DEBUG(Sc);
    }

    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPhaseChange::standardPhaseChange
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    phaseChangeModel(typeName, owner, dict),
    deltaMin_(readScalar(coeffDict_.lookup("deltaMin"))),
    L_(readScalar(coeffDict_.lookup("L"))),
    TbFactor_(coeffDict_.lookupOrDefault<scalar>("TbFactor", 1.1)),
    scaling_(coeffDict_.lookupOrDefault<scalar>("scaling",1.0)),
    dMassMax_(coeffDict_.lookupOrDefault<scalar>("dMassMax",5e-6))
{
    Info << "Mass transfer convective scaling set to " << scaling_ << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardPhaseChange::~standardPhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardPhaseChange::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy
)
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    // set local thermo properties
    const SLGThermo& thermo = film.thermo();
    const filmThermoModel& filmThermo = film.filmThermo();
    const label vapId = thermo.carrierId(filmThermo.name());

    // retrieve fields from film model
    const scalarField& delta = film.delta();
    const scalarField& YInf = film.YPrimary()[vapId];
    const scalarField& pInf = film.pPrimary();
    const scalarField& T = film.T();
    const scalarField& kappa = film.kappa();
    const scalarField& rho = film.rho();
    const scalarField& TInf = film.TPrimary();
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& muInf = film.muPrimary();
    const scalarField& magSf = film.magSf();
    const scalarField hInf(film.htcs().h());
    const vectorField dU(film.UPrimary() - film.Us());
    const filmRadiationModel& radiation = film.radiation();
    const scalarField Shs(radiation.ShsConst());
    // Info << "uprimary " << film.UPrimary();
    // Info << "us " << film.Us();
    const scalarField limMass
    (
        max(scalar(0.0), availableMass - deltaMin_*rho*magSf)
    );
    // const scalarField qRad(film.qRad());
    
    static scalarField dMassPrev(dMass.size(),0.0);

    scalar boilingMass = 0.0;
    scalar evaporationMass = 0.0;

    scalar ReMax = -GREAT;
    scalar ReMin =  GREAT;
    scalar ScMax = -GREAT;
    scalar ScMin =  GREAT;
    scalar dMassMax = -GREAT;
    scalar dMassMin =  GREAT;

    forAll(dMass, celli)
    {
        if (delta[celli] > deltaMin_)
        {
            // cell pressure [Pa]
            const scalar pc = max(min(pInf[celli],101325*1.1),101325*0.9);

            // calculate the boiling temperature
            const scalar Tb = filmThermo.Tb(pc);
            // const scalar Tb = 374.8;

            // estimate surface temperature based on energy balance
            const scalar d2k = delta[celli]/2.0/kappa[celli];
            const scalar Tf = T[celli];
            const scalar Tsurface = 
                (Tf+d2k*(hInf[celli]*TInf[celli]+Shs[celli]))
                /(1+d2k*hInf[celli]);

            // local temperature - impose lower limit of 200 K for stability
            // const scalar Tloc = min(TbFactor_*Tb, max(200.0, T[celli]));
            const scalar Tloc = min(TbFactor_*Tb, max(200.0, Tsurface));

            // saturation pressure [Pa]
            const scalar pSat = filmThermo.pv(pc, Tloc);

            // latent heat [J/kg]
            const scalar hVap = filmThermo.hl(pc, Tloc);

            // calculate mass transfer
            if (pSat >= 0.95*pc)
            {
                // boiling
                const scalar Cp = filmThermo.Cp(pc, Tloc);
                // const scalar Tcorr = 0.25*max(0.0, Tsurface - Tb);
                const scalar Tcorr = 0.25*max(0.0, Tloc - Tb);
                const scalar qCorr = limMass[celli]*Cp*(Tcorr);
                dMass[celli] = qCorr/hVap;
                boilingMass += dMass[celli];
            }
            else
            {
                // Primary region density [kg/m3]
                const scalar rhoInfc = rhoInf[celli];

                // Primary region viscosity [Pa.s]
                const scalar muInfc = muInf[celli];

                // Reynolds number
                scalar Re = rhoInfc*mag(dU[celli])*L_/muInfc;
                ReMax = max(ReMax,Re);
                ReMin = min(ReMin,Re);
                // limit Re to reasonable values, prevent runaway vaporization
                Re = min(1e6,Re);

                // molecular weight of vapour [kg/kmol]
                const scalar Wvap = thermo.carrier().W(vapId);

                // molecular weight of liquid [kg/kmol]
                const scalar Wliq = filmThermo.W();

                // vapour mass fraction at interface
                const scalar Ys = Wliq*pSat/(Wliq*pSat + Wvap*(pc - pSat));

                // vapour diffusivity [m2/s]
                const scalar Dab = filmThermo.D(pc, Tloc);

                // Schmidt number
                scalar Sc = muInfc/(rhoInfc*(Dab + ROOTVSMALL));
                if (Sc < 0.01)
                {
                    DEBUG(Sc);
                    DEBUG(muInfc);
                    DEBUG(rhoInfc);
                    DEBUG(Dab);
                    Sc = max(Sc,0.01);
                }
                ScMax = max(ScMax,Sc);
                ScMin = min(ScMin,Sc);

                // Sherwood number
                const scalar Sh = this->Sh(Re, Sc);

                // mass transfer coefficient [m/s]
                const scalar hm = scaling_*Sh*Dab/(L_ + ROOTVSMALL);
                // const scalar hm = scaling_;

                // add mass contribution to source
                dMass[celli] =
                    dt*magSf[celli]*rhoInfc*hm*(Ys - YInf[celli])/(1.0 - Ys);
                evaporationMass += dMass[celli];
                // Info << "Yinf " << YInf[celli] << endl;
                // Info << "Ys " << Ys << endl;
                // Info << "dMassb " << dMass[celli] << endl;
            }

            dMass[celli] = min(limMass[celli], max(0.0, dMass[celli]));
            // try to under-relax the vaporization rate
            dMass[celli] = min(dMassPrev[celli]*1.05+SMALL, max(0.0, dMass[celli]));
            // give a hard limit on vaporization rate, 
            // TODO: ideally should be based on cell volume and time step
            // dMass[celli] = min(5e-5, dMass[celli]);
            dMass[celli] = min(dMassMax_, dMass[celli]);
            // dMass[celli] = min(5e-7, dMass[celli]);
            dMassMax = max(dMassMax,dMass[celli]);
            dMassMin = min(dMassMin,dMass[celli]);

            // Info << "dMassa " << dMass[celli] << endl;
    
            dEnergy[celli] = dMass[celli]*hVap;
            dMassPrev[celli] = dMass[celli];
        }
    }
    reduce(ReMax,maxOp<scalar>());
    reduce(ReMin,minOp<scalar>());
    reduce(ScMax,maxOp<scalar>());
    reduce(ScMin,minOp<scalar>());
    reduce(dMassMax,maxOp<scalar>());
    // dMassMax_ = max(dMassMax,dMassMax_);
    reduce(dMassMin,minOp<scalar>());
    Info<<"standardPhaseChange::ReMax" << tab << film.time().timeName() << tab << ReMax <<endl;
    Info<<"standardPhaseChange::ReMin" << tab << film.time().timeName() << tab << ReMin <<endl;
    Info<<"standardPhaseChange::ScMax" << tab << film.time().timeName() << tab << ScMax <<endl;
    Info<<"standardPhaseChange::ScMin" << tab << film.time().timeName() << tab << ScMin <<endl;
    Info<<"standardPhaseChange::dMassMax" << tab << film.time().timeName() << tab << dMassMax <<endl;
    // Info<<"standardPhaseChange::dMassMax_" << tab << film.time().timeName() << tab << dMassMax_ <<endl;
    Info<<"standardPhaseChange::dMassMin" << tab << film.time().timeName() << tab << dMassMin <<endl;
    reduce(boilingMass,sumOp<scalar>());
    reduce(evaporationMass,sumOp<scalar>());
    Info << "boiling fraction " << boilingMass/(boilingMass+evaporationMass+SMALL) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

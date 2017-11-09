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

Class
    Foam::combustionModels::eddyDissipationFRVModel

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    eddyDissipationFRVModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationFRVModel_H
#define eddyDissipationFRVModel_H

#include "singleStepCombustion.H"

#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationFRVModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissipationFRVModel
:
    public singleStepCombustion<CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar C_;

        scalar Cd_;

        scalar Cstiff_;

        //- Sergey's Model

        //- Karlovitz Number
        volScalarField Ka_;
        volScalarField KaExt_;
        volScalarField KaMixed_;
        volScalarField KaExtMixed_;

        //- Flame Reactive Volume
        volScalarField FRV_;

        //- Ignition Factor
        volScalarField Fig_;

        //- Laminar Premixed Flame Speed (m/s)
        volScalarField SL_;
        volScalarField SLMixed_;

	//- Flame thickness
	volScalarField deltaFlame_;
	volScalarField deltaFlameMixed_;

	//- Reaction rates (W/m3)
	volScalarField dQFuel_;
	volScalarField dQFstar_;

	//- Heat Loss Term alpha
	volScalarField alphaLoss_;

	//- Expension Ratio
	volScalarField ExpR_;
	volScalarField ExpRMixed_;

	//- Gama = 6*G(n)*Kappa0*(Tb/T0)^1.75
	volScalarField Gama_;
	volScalarField GamaMixed_;

	//- Total Extinction Criteria
	volScalarField ExtNumber_;
	volScalarField ExtNumberMixed_;

	//- Laminar viscosity m^2/s
	volScalarField niuL_;

	//- Energy dissipation rate
	volScalarField epsSgs_;
	volScalarField epsG_;

	//- Adabatic Flame temperature
	volScalarField Tad_;
	volScalarField TadMixed_;

	//- Entrainment Ratio
	volScalarField ER_;

	//- Vairable Beta
	volScalarField Beta_;

	//- Reaction rates (kg/m3/s)
	volScalarField WFstar_;

	//- Scalars
	//- Time to Turnon the Extinction Model
	scalar tExt_;

	//- Zeldovich Number
	scalar ZN_;

	//- Flame velocity coefficient
	scalar TadAir_;
	scalar SLC1_;
	scalar SLC2_;

	//- Model calibration coefficients
	scalar cKa_;
	scalar cKapa_;
	scalar Cevap_;

	//- For adiabatic flame temperature model
	scalar XrExt_;
	scalar XO2Air_;
	scalar rN2_;
	scalar rCO2_;
	scalar rH2O_;
	scalar TFuel_;
	scalar TAir_;


    // Private Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;

        //- Disallow copy construct
        eddyDissipationFRVModel(const eddyDissipationFRVModel&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationFRVModel&);


public:

    //- Runtime type information
    TypeName("eddyDissipationFRVModel");


    // Constructors

        //- Construct from components
        eddyDissipationFRVModel
	(
	    const word& modelType, 
	    const fvMesh& mesh,
	    const word& combustionProperties,
	    const word& phaseName
	);


    //- Destructor
    virtual ~eddyDissipationFRVModel();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

	    //- Fuel consumption rate matrix                                     
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;                       

            //- Sensible enthalpy source term                                        
            //virtual tmp<volScalarField> Sh() const;
            //virtual tmp<volScalarField> ShFuel() const;
            //virtual tmp<volScalarField> ShFstar() const;
            virtual tmp<volScalarField> Qdot() const;
            virtual tmp<volScalarField> QdotFuel() const;
            virtual tmp<volScalarField> QdotFstar() const;

	    void calculateFlameTemperature();
	    void calculateReactiveVolume();


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyDissipationFRVModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

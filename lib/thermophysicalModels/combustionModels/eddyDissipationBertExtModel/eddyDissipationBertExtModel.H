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
    Foam::combustionModels::eddyDissipationBertExtModel

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    eddyDissipationBertExtModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationBertExtModel_H
#define eddyDissipationBertExtModel_H

#include "singleStepCombustion.H"

#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationBertExtModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissipationBertExtModel
:
    public singleStepCombustion<CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar C_;

        scalar Cd_;

        scalar Cstiff_;

        //- Bert's model

        //- Flame Extinction Factor
        volScalarField FEF_;

        //- Flame ReIgnition Factor
        volScalarField FIF_;

        //- Strain Rate (1/s)
        volScalarField strainRate_;

        //- Flame Extinction Temperature (K)
        volScalarField Textinction_;

        //- Adabatic Flame temperature
        volScalarField Tad_;

        //- Adabatic Flame temperature with heat loss
        volScalarField TadLoss_;

        //- Near wall cells
        volScalarField nearWallCells_;

        //- Reaction rates (kg/m3/s)
        volScalarField WCO2_;
        volScalarField WO2_;
        volScalarField WH2O_;
        volScalarField WF_;
        volScalarField WFstar_;
        volScalarField WFNet_;
        volScalarField WFstarNet_;

        scalar tExt_;
        scalar TFuelStarExt_;
        scalar TFuelExt_;
        scalar Cstrain_;
        scalar Cevap_;
        scalar XrExtinction_;
        bool nearWallExtinction_;
        bool radiativeHeatLoss_;


    // Private Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;

        //- Disallow copy construct
        eddyDissipationBertExtModel(const eddyDissipationBertExtModel&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationBertExtModel&);


public:

    //- Runtime type information
    TypeName("eddyDissipationBertExtModel");


    // Constructors

        //- Construct from components
        eddyDissipationBertExtModel
        (
            const word& modelType, 
            const fvMesh& mesh,
            const word& combustionProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~eddyDissipationBertExtModel();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix                                     
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;                       

            //- Sensible enthalpy source term                                        
            virtual tmp<volScalarField> Sh() const;

            void calculateFlameTemperature();


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyDissipationBertExtModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

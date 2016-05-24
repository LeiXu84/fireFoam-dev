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
    Foam::combustionModels::eddyDissMultiFuelStreamModel

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    eddyDissMultiFuelStreamModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissMultiFuelStreamModel_H
#define eddyDissMultiFuelStreamModel_H

#include "singleStepCombustion.H"

#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissMultiFuelStreamModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissMultiFuelStreamModel
:
    public singleStepCombustion<CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar C_;

        scalar Cd_;

        scalar Cstiff_;

        wordList fuelSpecies_;           // List of fuel species names in multi-fuel-streams model

        labelList fuelSpeciesID_;        // List of fuel species ID in the species table

        volScalarField fuelSpeciesSum_;  // Sum of mass-fractions of fuel species [ Can't use fuelID species to save this since that will violate sum(Yi) = 1]
                                         // This is used to distribite wFuel among fuel species appropriately..

    // Private Member Functions

        void initPars();    // Initializing the parameters 

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;

        //- Disallow copy construct
        eddyDissMultiFuelStreamModel(const eddyDissMultiFuelStreamModel&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissMultiFuelStreamModel&);


public:

    //- Runtime type information
    TypeName("eddyDissMultiFuelStreamModel");


    // Constructors

        //- Construct from components
        eddyDissMultiFuelStreamModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~eddyDissMultiFuelStreamModel();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //- Sensible enthalpy source term
            virtual tmp<volScalarField> Sh() const;


        // IO

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyDissMultiFuelStreamModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

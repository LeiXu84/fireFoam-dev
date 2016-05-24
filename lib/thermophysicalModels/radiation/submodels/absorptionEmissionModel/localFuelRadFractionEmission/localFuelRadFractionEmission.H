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

Class
    Foam::radiation::localFuelRadFractionEmission

Description
    This model computes emission as a fraction of the heat release for the continuous phase..
    The fraction is computed based on the relative local values of the mass-fractions for multiple streams of fuel (in the 
    context of the eddyDisMultiFuelStream combustion model).. 

    This is an optically thin model (i.e., no absorption is considered)..

SourceFiles
    localFuelRadFractionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef localFuelRadFractionEmission_H
#define localFuelRadFractionEmission_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "basicThermo.H"
#include "scalarIOList.H"
#include "wordIOList.H"

//#include "OpenFOAM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class localFuelRadFractionEmission Declaration
\*---------------------------------------------------------------------------*/

class localFuelRadFractionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Emission constant coefficient
        const scalar EhrrCoeff_;
        
        //- Switch for scaling radiant fraction
        Switch radScaling;
        
        //- Fuel species names for multiple streams of fuels..
        //wordList fuelStreams_;
        wordIOList fuelStreams_;

        //- Radiant fractions for different fuel streams
        //scalarList fuelRadFracs_; 
        scalarIOList fuelRadFracs_; 
        //scalarIOField fuelRadFracs_; 

        //- VolScalarField for radiant fraction
        volScalarField localRadFrac_;

public:

    //- Runtime type information
    TypeName("localFuelRadFractionEmission");


    // Constructors

        //- Construct from components
        localFuelRadFractionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~localFuelRadFractionEmission();


    // Member Operators

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::radiation::wsggmAbsorptionEmissionGreyCassol

Description
    Weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission,
    with Cassol et al coefficients.

    Grey formulation i.e. calculates total grey emissivity, grey absorption coefficient,
    to solve one RTE per ray.

    Look up dictionary file is in /constant. Example for use:

    /verbatim

    Pw_div_Pc   2; // use 0 for CO2 only, 1 for H2O only, 2 for Pw/Pc=2

    pathLength pathLength [ 0 1  0  0  0  0  0] 0.67;  // characteristic beam length in m
    
    Csoot                    1226; //- constant for soot absorption coefficient (coupling with volume fraction model)

    /end verbatim 
    
SourceFiles
    wsggmAbsorptionEmissionGreyCassol.C

\*---------------------------------------------------------------------------*/

#ifndef radiationwsggmAbsorptionEmissionGreyCassol_H
#define radiationwsggmAbsorptionEmissionGreyCassol_H

#include "absorptionEmissionModel.H"
//#include "hsReactionThermo.H"
//#include "hsCombustionThermo.H"
//#include "hReactionThermo.H"
#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H" 
#include "reactingMixture.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggmAbsorptionEmissionGreyCassol Declaration
\*---------------------------------------------------------------------------*/

class wsggmAbsorptionEmissionGreyCassol
:
    public absorptionEmissionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        volScalarField a_;

        //- Emission coefficient / [1/m]
        volScalarField e_;

        //- Emission contribution / [kg/(m s^3)]
        volScalarField E_;

        //characteristic path length / [m] from dict
        dimensionedScalar pathLength_;
        scalar Pw_div_Pc_;
        const scalar Csoot_;
        //- absorption coefficitent for individual grey gases
        scalarListList emissivityCoeffs_low;
        scalarListList emissivityCoeffs_mid;
        scalarListList emissivityCoeffs_high;
        //- fitting Factors
        List<scalarListList> fittingFactors_low;
        List<scalarListList> fittingFactors_mid;
        List<scalarListList> fittingFactors_high;


        
        // absorption/emission coeff


public:

    //- Runtime type information
    TypeName("wsggmAbsorptionEmissionGreyCassol");


    // Constructors

        //- Construct from components
        wsggmAbsorptionEmissionGreyCassol
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Destructor
    virtual ~wsggmAbsorptionEmissionGreyCassol();


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
            return false; //true (modified by Ivan Sikic 01/10/2014)
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

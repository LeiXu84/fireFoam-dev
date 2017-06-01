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

\*---------------------------------------------------------------------------*/

#include "wsggmAbsorptionEmissionGreyCassol.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wsggmAbsorptionEmissionGreyCassol, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wsggmAbsorptionEmissionGreyCassol,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmissionGreyCassol::wsggmAbsorptionEmissionGreyCassol
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("e", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
//    emissivityCoeffs_(coeffsDict_.lookup("emissivityCoeffs")),
//    fittingFactors_(coeffsDict_.lookup("fittingFactors")),
    pathLength_(coeffsDict_.lookup("pathLength")),
//  Pw_div_Pc_ =2, 1, 1000(Pc->0atm)
    Pw_div_Pc_(readScalar(coeffsDict_.lookup("Pw_div_Pc"))),
    Csoot_(readScalar(coeffsDict_.lookup("Csoot")))
{

//  low temperature 400-2500; NO DATA FOR HIGHER TEMPS.

	  emissivityCoeffs_low.setSize(3);
	  forAll(emissivityCoeffs_low,i)
	     {
	      emissivityCoeffs_low[i].setSize(4);
	     }
// first dimension: 0 for Pw/Pc=2; 1 for H2O only; 2 for CO2 only
          emissivityCoeffs_low[0][0]=0.192;
          emissivityCoeffs_low[0][1]=1.719;
          emissivityCoeffs_low[0][2]=11.37;
          emissivityCoeffs_low[0][3]=111.016;

          emissivityCoeffs_low[1][0]=0.171;
          emissivityCoeffs_low[1][1]=1.551;
          emissivityCoeffs_low[1][2]=5.562;
          emissivityCoeffs_low[1][3]=49.159;

          emissivityCoeffs_low[2][0]=0.138;
          emissivityCoeffs_low[2][1]=1.895;
          emissivityCoeffs_low[2][2]=13.301;
          emissivityCoeffs_low[2][3]=340.811;

	  fittingFactors_low.setSize(3);
	  forAll(fittingFactors_low,i)
	    {
	    	fittingFactors_low[i].setSize(4);
	    	 for(label j=0; j<4;j++)
	    	  {
	    	  	fittingFactors_low[i][j].setSize(5);
	    	  }
	    }

// Pw/Pc = 2
	  fittingFactors_low[0][0][0]=0.05617;
	  fittingFactors_low[0][0][1]=78.44;
	  fittingFactors_low[0][0][2]=-85.63;
	  fittingFactors_low[0][0][3]=42.46;
	  fittingFactors_low[0][0][4]=-74.4;

	  fittingFactors_low[0][1][0]=0.1426;
	  fittingFactors_low[0][1][1]=17.95;
	  fittingFactors_low[0][1][2]=-1.077;
	  fittingFactors_low[0][1][3]=-6.971;
	  fittingFactors_low[0][1][4]=17.74;

	  fittingFactors_low[0][2][0]=0.1362;
	  fittingFactors_low[0][2][1]=25.74;
	  fittingFactors_low[0][2][2]=-37.11;
	  fittingFactors_low[0][2][3]=15.7;
	  fittingFactors_low[0][2][4]=-22.67;

	  fittingFactors_low[0][3][0]=0.1222;
	  fittingFactors_low[0][3][1]=-2.327;
	  fittingFactors_low[0][3][2]=-7.492;
	  fittingFactors_low[0][3][3]=4.275;
	  fittingFactors_low[0][3][4]=-6.608;
// H2O only
	  fittingFactors_low[1][0][0]=0.06617;
	  fittingFactors_low[1][0][1]=55.48;
	  fittingFactors_low[1][0][2]=-48.41;
	  fittingFactors_low[1][0][3]=22.27;
	  fittingFactors_low[1][0][4]=-40.17;

	  fittingFactors_low[1][1][0]=0.11045;
	  fittingFactors_low[1][1][1]=0.576;
	  fittingFactors_low[1][1][2]=24;
	  fittingFactors_low[1][1][3]=-17.01;
	  fittingFactors_low[1][1][4]=30.96;

	  fittingFactors_low[1][2][0]=-0.04915;
	  fittingFactors_low[1][2][1]=70.63;
	  fittingFactors_low[1][2][2]=-70.12;
	  fittingFactors_low[1][2][3]=26.07;
	  fittingFactors_low[1][2][4]=-34.94;

	  fittingFactors_low[1][3][0]=0.23675;
	  fittingFactors_low[1][3][1]=-18.91;
	  fittingFactors_low[1][3][2]=-0.907;
	  fittingFactors_low[1][3][3]=4.082;
	  fittingFactors_low[1][3][4]=-8.778;
// CO2 only
	  fittingFactors_low[2][0][0]=0.0999;
	  fittingFactors_low[2][0][1]=64.41;
	  fittingFactors_low[2][0][2]=-86.94;
	  fittingFactors_low[2][0][3]=41.27;
	  fittingFactors_low[2][0][4]=-67.74;

	  fittingFactors_low[2][1][0]=0.00942;
	  fittingFactors_low[2][1][1]=10.36;
	  fittingFactors_low[2][1][2]=-2.277;
	  fittingFactors_low[2][1][3]=-2.134;
	  fittingFactors_low[2][1][4]=6.497;

	  fittingFactors_low[2][2][0]=0.14511;
	  fittingFactors_low[2][2][1]=-30.73;
	  fittingFactors_low[2][2][2]=37.65;
	  fittingFactors_low[2][2][3]=-18.41;
	  fittingFactors_low[2][2][4]=30.16;

	  fittingFactors_low[2][3][0]=-0.02915;
	  fittingFactors_low[2][3][1]=25.23;
	  fittingFactors_low[2][3][2]=-26.1;
	  fittingFactors_low[2][3][3]=9.965;
	  fittingFactors_low[2][3][4]=-13.26;


	  forAll(fittingFactors_low,i)
	     {
	     	 for(label j=0; j<4;j++)
	     	  {
	     	  	  fittingFactors_low[i][j][0] /=1;
	     	  	  fittingFactors_low[i][j][1] /=1.0e5;
	     	  	  fittingFactors_low[i][j][2] /=1.0e8;
	     	  	  fittingFactors_low[i][j][3] /=1.0e11;
	     	  	  fittingFactors_low[i][j][4] /=1.0e15;
	     	  }
	     }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmissionGreyCassol::~wsggmAbsorptionEmissionGreyCassol()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmissionGreyCassol::aCont(const label bandI) const
{

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );

    const psiReactionThermo& thermo= mesh_.lookupObject<psiReactionThermo>("thermophysicalProperties");

    //access species mass fractions
    const PtrList<volScalarField>& Y = thermo.composition().Y();
    //access specie thermo data
    const PtrList<gasHThermoPhysics> & specieThermo =
        dynamic_cast<const reactingMixture<gasHThermoPhysics>&>  (thermo).speciesData();
    // get index of CO2 in mixture
    const volScalarField& fv =mesh_.lookupObject<volScalarField>("fv");
    label indexCO2= dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["CO2"];
    // get index of H2O in mixture
    label indexH2O= dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["H2O"];

    volScalarField emissivity(
            IOobject
            (
                "emissivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
    volScalarField pressurePathLength(
            IOobject
            (
                "pressurePathLength",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
    volScalarField weightingFactor(
            IOobject
            (
                "weightingFactor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
//    const scalar wCO2=44;
//    const scalar wH2O=18;

    volScalarField wMean(
            IOobject
            (
                "wMean",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0) // kg/kmol
        );

    forAll(Y,specieI)
    {
        wMean+=Y[specieI]/specieThermo[specieI].W();
    }

    wMean=1/wMean;

    pressurePathLength=wMean*(thermo.p()/101325/dimensionedScalar("unity",dimPressure,1.0))*(Y[indexCO2]/specieThermo[indexCO2].W()
                      +Y[indexH2O]/specieThermo[indexH2O].W())*pathLength_.value();

    volScalarField limitedDimlessTemperature = min(thermo.T()/dimensionedScalar("unity",dimTemperature, 1.0), 3000.0);




    if(Pw_div_Pc_==2) // Mixture : Pw/Pc =2
    {
       forAll(limitedDimlessTemperature, cellI)
        {
    	  if(limitedDimlessTemperature[cellI]<=2500)
    	    {
    	      forAll(emissivityCoeffs_low[0],i)
               {
                  weightingFactor[cellI] = 0.0;
                  for(label j=0; j<fittingFactors_low[0][0].size();j++)
                    {
                      weightingFactor[cellI] += fittingFactors_low[0][i][j]*pow(limitedDimlessTemperature[cellI], j);
                    }
                  emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_low[0][i] *pressurePathLength[cellI]));
               }
            }

        }
     }

    if(Pw_div_Pc_==1) // H2O ONLY
    {
       forAll(limitedDimlessTemperature, cellI)
        {
    	  if(limitedDimlessTemperature[cellI]<=2500)
    	     {
    	       forAll(emissivityCoeffs_low[1],i)
                 {
                    weightingFactor[cellI] = 0.0;
                    for(label j=0; j<fittingFactors_low[1][0].size();j++)
                     {
                       weightingFactor[cellI] += fittingFactors_low[1][i][j]*pow(limitedDimlessTemperature[cellI], j);
                     }
                    emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_low[1][i] * pressurePathLength[cellI]));
                 }
             }

        }
     }

    if(Pw_div_Pc_==0) // CO2 ONLY
    {
      forAll(limitedDimlessTemperature, cellI)
    	{
          if(limitedDimlessTemperature[cellI]<=2500)
    	    {
    	      forAll(emissivityCoeffs_low[2],i)
                {
                   weightingFactor[cellI] = 0.0;
                   for(label j=0; j<fittingFactors_low[2][0].size();j++)
                    {
                      weightingFactor[cellI] += fittingFactors_low[2][i][j]*pow(limitedDimlessTemperature[cellI], j);
                    }
                    emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_low[2][i] * pressurePathLength[cellI]));
                }
            }

        }
    }

   volScalarField aSoot = Csoot_*fv*limitedDimlessTemperature;
   aSoot.dimensions().reset(dimensionSet(0,-1,0,0,0,0,0));
   // Info << "weightingFactor = "  << weightingFactor << endl;


   emissivity = min(emissivity,0.9999);
//    a_= (-1)* log(1-emissivity) / pathLength_;

    ta.ref()= (-1)* log(1-emissivity) / pathLength_ + aSoot;

    return ta;

}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmissionGreyCassol::eCont(const label bandI) const
{

    return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmissionGreyCassol::ECont(const label bandI) const
{
    tmp<volScalarField> tE=E_;

    return tE;

}

// ************************************************************************* //

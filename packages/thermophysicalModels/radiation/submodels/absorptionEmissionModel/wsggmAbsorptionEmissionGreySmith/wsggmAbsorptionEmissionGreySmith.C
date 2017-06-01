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

#include "wsggmAbsorptionEmissionGreySmith.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wsggmAbsorptionEmissionGreySmith, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wsggmAbsorptionEmissionGreySmith,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmissionGreySmith::wsggmAbsorptionEmissionGreySmith
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


//  low temperature 600-2400; middle temperature 2400-2500; high temperature 2500-3000.
	  emissivityCoeffs_low.setSize(5); // modified by Ivan Sikic 13/08/2014
	  forAll(emissivityCoeffs_low,i)
	     {
	      emissivityCoeffs_low[i].setSize(3);
	     }
// first dimension: 0 for Pw/Pc=2; 1 for Pw/Pc=1; 2 for Pc->0atm
          emissivityCoeffs_low[0][0]=0.4201;
          emissivityCoeffs_low[0][1]=6.516;
          emissivityCoeffs_low[0][2]=131.9;
          emissivityCoeffs_low[1][0]=0.4303;
          emissivityCoeffs_low[1][1]=7.055;
          emissivityCoeffs_low[1][2]=178.1;
          emissivityCoeffs_low[2][0]=0.3966;
          emissivityCoeffs_low[2][1]=15.64;
          emissivityCoeffs_low[2][2]=394.3;
// modified by Ivan Sikic 14/08/2014: 3 for Pw->0atm (Pw_div_Pc=0.001), 4 for Pw=1atm (Pw_div_Pc=10)
	  emissivityCoeffs_low[3][0]=0.4098;
	  emissivityCoeffs_low[3][1]=6.325;
	  emissivityCoeffs_low[3][2]=120.5;
	  emissivityCoeffs_low[4][0]=0.4496;
	  emissivityCoeffs_low[4][1]=7.113;
	  emissivityCoeffs_low[4][2]=119.7;
// end of modification

	  fittingFactors_low.setSize(5); // size was 3
	  forAll(fittingFactors_low,i)
	    {
	    	fittingFactors_low[i].setSize(3);
	    	 for(label j=0; j<3;j++)
	    	  {
	    	  	fittingFactors_low[i][j].setSize(4);
	    	  }
	    }


	  fittingFactors_low[0][0][0]=6.508;
	  fittingFactors_low[0][0][1]=-5.551;
	  fittingFactors_low[0][0][2]=3.029;
	  fittingFactors_low[0][0][3]=-5.353;
	  fittingFactors_low[0][1][0]=-0.2504;
	  fittingFactors_low[0][1][1]=6.112;
	  fittingFactors_low[0][1][2]=-3.882;
	  fittingFactors_low[0][1][3]=6.528;
	  fittingFactors_low[0][2][0]=2.718;
	  fittingFactors_low[0][2][1]=-3.118;
	  fittingFactors_low[0][2][2]=1.221;
	  fittingFactors_low[0][2][3]=-1.612;

	  fittingFactors_low[1][0][0]=5.150;
	  fittingFactors_low[1][0][1]=-2.303;
	  fittingFactors_low[1][0][2]=0.9779;
	  fittingFactors_low[1][0][3]=-1.494;
	  fittingFactors_low[1][1][0]=0.7749;
	  fittingFactors_low[1][1][1]=3.399;
	  fittingFactors_low[1][1][2]=-2.297;
	  fittingFactors_low[1][1][3]=3.770;
	  fittingFactors_low[1][2][0]=1.907;
	  fittingFactors_low[1][2][1]=-1.824;
	  fittingFactors_low[1][2][2]=0.5608;
	  fittingFactors_low[1][2][3]=-0.5122;

	  fittingFactors_low[2][0][0]=0.4334;
	  fittingFactors_low[2][0][1]=2.620;
	  fittingFactors_low[2][0][2]=-1.560;
	  fittingFactors_low[2][0][3]=2.565;
	  fittingFactors_low[2][1][0]=-0.4814;
	  fittingFactors_low[2][1][1]=2.822;
	  fittingFactors_low[2][1][2]=-1.794;
	  fittingFactors_low[2][1][3]=3.274;
	  fittingFactors_low[2][2][0]=0.5492;
	  fittingFactors_low[2][2][1]=0.1087;
	  fittingFactors_low[2][2][2]=-0.3500;
	  fittingFactors_low[2][2][3]=0.9123;
// modified by Ivan Sikic 14/08/2014
	  fittingFactors_low[3][0][0]=5.977;
	  fittingFactors_low[3][0][1]=-5.119;
	  fittingFactors_low[3][0][2]=3.042;
	  fittingFactors_low[3][0][3]=-5.564;
	  fittingFactors_low[3][1][0]=0.5677;
	  fittingFactors_low[3][1][1]=3.333;
	  fittingFactors_low[3][1][2]=-1.967;
	  fittingFactors_low[3][1][3]=2.718;
	  fittingFactors_low[3][2][0]=1.8;
	  fittingFactors_low[3][2][1]=-2.334;
	  fittingFactors_low[3][2][2]=1.008;
	  fittingFactors_low[3][2][3]=-1.454;

	  fittingFactors_low[4][0][0]=6.324;
	  fittingFactors_low[4][0][1]=-8.358;
	  fittingFactors_low[4][0][2]=6.135;
	  fittingFactors_low[4][0][3]=-13.03;
	  fittingFactors_low[4][1][0]=-0.2016;
	  fittingFactors_low[4][1][1]=7.145;
	  fittingFactors_low[4][1][2]=-5.212;
	  fittingFactors_low[4][1][3]=9.868;
	  fittingFactors_low[4][2][0]=3.5;
	  fittingFactors_low[4][2][1]=-5.04;
	  fittingFactors_low[4][2][2]=2.425;
	  fittingFactors_low[4][2][3]=-3.888;
// end of modification

	  forAll(fittingFactors_low,i)
	     {
	     	 for(label j=0; j<3;j++)
	     	  {
	     	  	  fittingFactors_low[i][j][0] /=10;
	     	  	  fittingFactors_low[i][j][1] /=1.0e4;
	     	  	  fittingFactors_low[i][j][2] /=1.0e7;
	     	  	  fittingFactors_low[i][j][3] /=1.0e11;
	     	  }
	     }
// 2400-2500
	  emissivityCoeffs_mid.setSize(3);
	  forAll(emissivityCoeffs_mid,i)
	     {
	      emissivityCoeffs_mid[i].setSize(3);
	     }
          emissivityCoeffs_mid[0][0]=0.527;
          emissivityCoeffs_mid[0][1]=3.78;
          emissivityCoeffs_mid[0][2]=99.54;
          emissivityCoeffs_mid[1][0]=0.464;
          emissivityCoeffs_mid[1][1]=3.47;
          emissivityCoeffs_mid[1][2]=121.6;
          emissivityCoeffs_mid[2][0]=0.3966;
          emissivityCoeffs_mid[2][1]=15.64;
          emissivityCoeffs_mid[2][2]=394.3;
	  fittingFactors_mid.setSize(3);
	  forAll(fittingFactors_mid,i)
	    {
	    	fittingFactors_mid[i].setSize(3);
	    	 for(label j=0; j<3;j++)
	    	  {
	    	  	fittingFactors_mid[i][j].setSize(4);
	    	  }
	    }

	  fittingFactors_mid[0][0][0]=0.132;
	  fittingFactors_mid[0][0][1]=0.0000725;
	  fittingFactors_mid[0][0][2]=0.0;
	  fittingFactors_mid[0][0][3]=0.0;
	  fittingFactors_mid[0][1][0]=0.547;
	  fittingFactors_mid[0][1][1]=-0.000171;
	  fittingFactors_mid[0][1][2]=-0.0;
	  fittingFactors_mid[0][1][3]=-0.0;
	  fittingFactors_mid[0][2][0]=0.0489;
	  fittingFactors_mid[0][2][1]=-0.0000176;
	  fittingFactors_mid[0][2][2]=0.0;
	  fittingFactors_mid[0][2][3]=0.0;

	  fittingFactors_mid[1][0][0]=0.136;
	  fittingFactors_mid[1][0][1]=0.0000726;
	  fittingFactors_mid[1][0][2]=0.0;
	  fittingFactors_mid[1][0][3]=0.0;
	  fittingFactors_mid[1][1][0]=0.516;
	  fittingFactors_mid[1][1][1]=-0.000163;
	  fittingFactors_mid[1][1][2]=0.0;
	  fittingFactors_mid[1][1][3]=0.0;
	  fittingFactors_mid[1][2][0]=0.0517;
	  fittingFactors_mid[1][2][1]=-0.0000176;
	  fittingFactors_mid[1][2][2]=0.0;
	  fittingFactors_mid[1][2][3]=0.0;

	  fittingFactors_mid[2][0][0]=0.4334;
	  fittingFactors_mid[2][0][1]=2.620;
	  fittingFactors_mid[2][0][2]=-1.560;
	  fittingFactors_mid[2][0][3]=2.565;
	  fittingFactors_mid[2][1][0]=-0.4814;
	  fittingFactors_mid[2][1][1]=2.822;
	  fittingFactors_mid[2][1][2]=-1.794;
	  fittingFactors_mid[2][1][3]=3.274;
	  fittingFactors_mid[2][2][0]=0.5492;
	  fittingFactors_mid[2][2][1]=0.1087;
	  fittingFactors_mid[2][2][2]=-0.3500;
	  fittingFactors_mid[2][2][3]=0.9123;
//	  i=2;
	  for(label j=0; j<3;j++)
	    {
	    	  fittingFactors_mid[2][j][0] /=10;
	    	  fittingFactors_mid[2][j][1] /=1.0e4;
	    	  fittingFactors_mid[2][j][2] /=1.0e7;
	    	  fittingFactors_mid[2][j][3] /=1.0e11;
	    }
//  2500-3000
	  emissivityCoeffs_high.setSize(3);
	  forAll(emissivityCoeffs_high,i)
	     {
	      emissivityCoeffs_high[i].setSize(3);
	     }
          emissivityCoeffs_high[0][0]=0.527;
          emissivityCoeffs_high[0][1]=3.78;
          emissivityCoeffs_high[0][2]=99.54;
          emissivityCoeffs_high[1][0]=0.464;
          emissivityCoeffs_high[1][1]=3.47;
          emissivityCoeffs_high[1][2]=121.6;
          emissivityCoeffs_high[2][0]=0.3966;
          emissivityCoeffs_high[2][1]=15.64;
          emissivityCoeffs_high[2][2]=394.3;
	  fittingFactors_high.setSize(3);
	  forAll(fittingFactors_high,i)
	    {
	    	fittingFactors_high[i].setSize(3);
	    	 for(label j=0; j<3;j++)
	    	  {
	    	  	fittingFactors_high[i][j].setSize(4);
	    	  }
	    }

	  fittingFactors_high[0][0][0]=0.430;
	  fittingFactors_high[0][0][1]=-0.0000472;
	  fittingFactors_high[0][0][2]=0.0;
	  fittingFactors_high[0][0][3]=0.0;
	  fittingFactors_high[0][1][0]=0.37;
	  fittingFactors_high[0][1][1]=-0.000101;
	  fittingFactors_high[0][1][2]=0.0;
	  fittingFactors_high[0][1][3]=0.0;
	  fittingFactors_high[0][2][0]=0.0184;
	  fittingFactors_high[0][2][1]=-0.00000511;
	  fittingFactors_high[0][2][2]=0.0;
	  fittingFactors_high[0][2][3]=0.0;

	  fittingFactors_high[1][0][0]=0.464;
	  fittingFactors_high[1][0][1]=-0.0000596;
	  fittingFactors_high[1][0][2]=0.0;
	  fittingFactors_high[1][0][3]=0.0;
	  fittingFactors_high[1][1][0]=0.336;
	  fittingFactors_high[1][1][1]=-0.0000909;
	  fittingFactors_high[1][1][2]=0.0;
	  fittingFactors_high[1][1][3]=0.0;
	  fittingFactors_high[1][2][0]=0.0245;
	  fittingFactors_high[1][2][1]=-0.00000654;
	  fittingFactors_high[1][2][2]=0.0;
	  fittingFactors_high[1][2][3]=0.0;

	  fittingFactors_high[2][0][0]=0.4334;
	  fittingFactors_high[2][0][1]=2.620;
	  fittingFactors_high[2][0][2]=-1.560;
	  fittingFactors_high[2][0][3]=2.565;
	  fittingFactors_high[2][1][0]=-0.4814;
	  fittingFactors_high[2][1][1]=2.822;
	  fittingFactors_high[2][1][2]=-1.794;
	  fittingFactors_high[2][1][3]=3.274;
	  fittingFactors_high[2][2][0]=0.5492;
	  fittingFactors_high[2][2][1]=0.1087;
	  fittingFactors_high[2][2][2]=-0.3500;
	  fittingFactors_high[2][2][3]=0.9123;

	  for(label j=0; j<3;j++)
	    {
	    	  fittingFactors_high[2][j][0] /=10;
	    	  fittingFactors_high[2][j][1] /=1.0e4;
	    	  fittingFactors_high[2][j][2] /=1.0e7;
	    	  fittingFactors_high[2][j][3] /=1.0e11;
	    }


	}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmissionGreySmith::~wsggmAbsorptionEmissionGreySmith()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmissionGreySmith::aCont(const label bandI) const
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
    const volScalarField& fv =mesh_.lookupObject<volScalarField>("fv"); // get calculated soot vol fraction

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




    if(Pw_div_Pc_==2)
    {
       forAll(limitedDimlessTemperature, cellI)
        {
    	  if(limitedDimlessTemperature[cellI]<=2400)
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
    	  else if(limitedDimlessTemperature[cellI]<=2500)
            {
              forAll(emissivityCoeffs_mid[0],i)
               {
                  weightingFactor[cellI] = 0.0;
                  for(label j=0; j<fittingFactors_mid[0][0].size();j++)
                    {
                      weightingFactor[cellI] += fittingFactors_mid[0][i][j]*pow(limitedDimlessTemperature[cellI], j);
                    }
                  emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_mid[0][i] *  pressurePathLength[cellI]));
               }
            }
    	  else if(limitedDimlessTemperature[cellI]<=3000)
    	    {
              forAll(emissivityCoeffs_high[0],i)
               {
                  weightingFactor[cellI] = 0.0;
                  for(label j=0; j<fittingFactors_high[0][0].size();j++)
                    {
                      weightingFactor[cellI] += fittingFactors_high[0][i][j]*pow(limitedDimlessTemperature[cellI], j);
                    }
                  emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_high[0][i] * pressurePathLength[cellI]));
               }
            }
        }
     }

    if(Pw_div_Pc_==1)
    {
       forAll(limitedDimlessTemperature, cellI)
        {
    	  if(limitedDimlessTemperature[cellI]<=2400)
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
    	   else if(limitedDimlessTemperature[cellI]<=2500)
    	     {
    	       forAll(emissivityCoeffs_mid[1],i)
                 {
                    weightingFactor[cellI] = 0.0;
                    for(label j=0; j<fittingFactors_mid[1][0].size();j++)
                     {
                       weightingFactor[cellI] += fittingFactors_mid[1][i][j]*pow(limitedDimlessTemperature[cellI], j);
                     }
                    emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_mid[1][i] * pressurePathLength[cellI]));
                 }
             }
           else if(limitedDimlessTemperature[cellI]<=3000)
    	     {
    	       forAll(emissivityCoeffs_high[1],i)
                 {
                    weightingFactor[cellI] = 0.0;
                    for(label j=0; j<fittingFactors_high[1][0].size();j++)
                     {
                       weightingFactor[cellI] += fittingFactors_high[1][i][j]*pow(limitedDimlessTemperature[cellI], j);
                     }
                    emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_high[1][i] * pressurePathLength[cellI]));
                 }
             }
        }
     }

    if(Pw_div_Pc_==0.001) // CO2 only, Pc -> 0
    {
      forAll(limitedDimlessTemperature, cellI)
    	{
          if(limitedDimlessTemperature[cellI]<=2400)
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
    	  else if(limitedDimlessTemperature[cellI]<=2500)
    	    {
    	      forAll(emissivityCoeffs_mid[2],i)
                {
                   weightingFactor[cellI] = 0.0;
                   for(label j=0; j<fittingFactors_mid[2][0].size();j++)
                    {
                      weightingFactor[cellI] += fittingFactors_mid[2][i][j]*pow(limitedDimlessTemperature[cellI], j);
                    }
                   emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_mid[2][i] * pressurePathLength[cellI]));
                }
            }
    	 else if(limitedDimlessTemperature[cellI]<=3000)
    	   {
    	     forAll(emissivityCoeffs_high[2],i)
               {
                   weightingFactor[cellI] = 0.0;
                   for(label j=0; j<fittingFactors_high[2][0].size();j++)
                    {
                      weightingFactor[cellI] += fittingFactors_high[2][i][j]*pow(limitedDimlessTemperature[cellI], j);
                    }
                   emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_high[2][i] * pressurePathLength[cellI]));
               }
           }
        }
    }

// modified by Ivan Sikic 14/08/2014

    if(Pw_div_Pc_==1000) // H2O only, Pw -> 0
    {
       forAll(limitedDimlessTemperature, cellI)
        {
    	  if(limitedDimlessTemperature[cellI]<=2400)
    	     {
    	       forAll(emissivityCoeffs_low[3],i)
                 {
                    weightingFactor[cellI] = 0.0;
                    for(label j=0; j<fittingFactors_low[3][0].size();j++)
                     {
                       weightingFactor[cellI] += fittingFactors_low[3][i][j]*pow(limitedDimlessTemperature[cellI], j);
                     }
                    emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_low[3][i] * pressurePathLength[cellI]));

                  }
             }
        }
    }

    if(Pw_div_Pc_==10) // H2O only, Pw = 1atm
    {
       forAll(limitedDimlessTemperature, cellI)
        {
    	  if(limitedDimlessTemperature[cellI]<=2400)
    	     {
    	       forAll(emissivityCoeffs_low[4],i)
                 {
                    weightingFactor[cellI] = 0.0;
                    for(label j=0; j<fittingFactors_low[4][0].size();j++)
                     {
                       weightingFactor[cellI] += fittingFactors_low[4][i][j]*pow(limitedDimlessTemperature[cellI], j);
                     }
                    emissivity[cellI] += weightingFactor[cellI] * (1 - exp( (-1) *emissivityCoeffs_low[4][i] * pressurePathLength[cellI]));
                 }
             }
        }
    }

   volScalarField aSoot = Csoot_*fv*limitedDimlessTemperature;
   aSoot.dimensions().reset(dimensionSet(0,-1,0,0,0,0,0));

   emissivity = min(emissivity,0.9999);
//    a_= (-1)* log(1-emissivity) / pathLength_;

    ta.ref() = (-1)* log(1-emissivity) / pathLength_ + aSoot;

    return ta;

}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmissionGreySmith::eCont(const label bandI) const
{

    return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmissionGreySmith::ECont(const label bandI) const
{
    tmp<volScalarField> tE=E_;

    return tE;

}

// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

#include "spSootRadFracEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "zeroGradientFvPatchFields.H"
#include "basicMultiComponentMixture.H"
#include "turbulentFluidThermoModel.H"

#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(spSootRadFracEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            spSootRadFracEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::spSootRadFracEmission::spSootRadFracEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    thermo_(mesh.lookupObject<basicThermo>("thermophysicalProperties")),
    turbulence_(mesh.lookupObject<turbulenceModel>("turbulenceProperties")),
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
    radScaling(coeffsDict_.lookupOrDefault<Switch>("radScaling",false)),
    Ehrr1_(coeffsDict_.lookupOrDefault<scalar>("Ehrr1",0.3)),
    Ehrr2_(coeffsDict_.lookupOrDefault<scalar>("Ehrr2",0.3)),
    patchName1_(coeffsDict_.lookup("patch1")),
    patchName2_(coeffsDict_.lookup("patch2")),
    YO2Inf(coeffsDict_.lookupOrDefault<scalar>("YO2Inf",0.23301)),
    Ceta(coeffsDict_.lookupOrDefault<scalar>("Ceta",0.04)),
    Ceta0(coeffsDict_.lookupOrDefault<scalar>("Ceta0",4)), 
    tableName_("none"),
    globalYO2_(coeffsDict_.lookupOrDefault<bool>("globalYO2","true")),
    fuel_("none"),
    tableOxyMassFracs_(),
    psiRTables_(),
    sr_(),
    hp_(),
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
        dimensionedScalar("zero", pow(dimTime,-1), 1e3)
    )



{

    readTableData();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::spSootRadFracEmission::~spSootRadFracEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::spSootRadFracEmission::aCont(const label bandI) const
{

    tmp<volScalarField> a
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    return a;

}


Foam::tmp<Foam::volScalarField>
Foam::radiation::spSootRadFracEmission::eCont(const label bandI) const
{
    tmp<volScalarField> e
    (
        new volScalarField
        (
            IOobject
            (
                "eCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    return e;
}


void Foam::radiation::spSootRadFracEmission::readTableData()
{
    Info << "Reading soot-radiation lookup table data\n";

    // Provide tableName in radiationPropertiesDict spSootRadFracEmissionCoeffs section
    // the tableName will correspond to the file where lookup table values are stored
    tableName_ = this->coeffsDict_.subDict("lookupTableCoeffs").template lookupOrDefault< word >("tableName","");
    fileName constant;
    if(Pstream::parRun()){
        constant =
            ".."/mesh().time().constant();
    }
    else{
        //Info << "tableName_: " << tableName_ << endl;
        constant =
            mesh().time().constant();
    }

    {
        IOdictionary dict
        (
            IOobject
            (
                tableName_,
                constant,
                mesh(), 
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

	fuel_ = dict.template lookupOrDefault< word >("fuel","not defined");
	//Info << "Fuel " << fuel_<< nl;
	dict.lookup("SR") >> sr_;
	dict.lookup("Hp") >> hp_;
    	scalar size = sr_.size()*hp_.size();
	//Info << "sr_" << sr_ << nl;
	//Info << "hp_" << hp_ << nl;

	const dictionary& oxyMassFracsDict = dict.subDict("oxyMassFracs");
	//Info << oxyMassFracsDict.toc() << nl;
	const wordList oxyMassFracNames(oxyMassFracsDict.toc());
	//Info << "oxyMassFracNames = " << oxyMassFracNames << nl;
	forAll(oxyMassFracNames,i)
	{
		const dictionary& oxyMassFracSubDict = oxyMassFracsDict.subDict(oxyMassFracNames[i]);
		tableOxyMassFracs_.append(readScalar(oxyMassFracSubDict.lookup("oxyMassFrac")));
		//Info << "oxy mass frac\n" << nl;
		//Info << "tableOxyMassFracs_: " << tableOxyMassFracs_ << nl;
		//Info << "tableOxyMassFracs_.size(): " << tableOxyMassFracs_.size() << nl;

		psiRTables_.append(oxyMassFracSubDict.lookup("PsiR"));
		//Info << "psiRTables[" << i << "]: " << psiRTables_[i] << nl;
		//Info << "psiRTables[" << i << "].size(): " << psiRTables_[i].size() << nl;
            	if(psiRTables_[i].size() != size)
            	{
            	    FatalErrorIn("spSootRadFracEmission::readTableData()")
            	        << "sr_.size()*hp_.size() not equal to psiRTables[" << i << "].size()" << nl
            	        << abort(FatalError);
            	}
	}

	
	//sr_.resize(srTables_[0].size());
	//hp_.resize(hpTables_[0].size());
	//psiR_.resize(psiRTables_[0].size());
	//sr_ = srTables_[0];
	//hp_ = hpTables_[0];
	//psiR_ = psiRTables_[0];
	//Info << "SR = " << sr_ << nl;
	//Info << "Hp = " << hp_ << nl;
	//Info << "PsiR = " << psiR_ << nl;
    }

    // Computing Ceta as a function of Ck
    // Ceta = Ceta0 x 10^-3 / Ck
    const compressible::LESModel& lesModel = mesh_.lookupObject<compressible::LESModel>(turbulenceModel::propertiesName); //("LESProperties");
    scalar Ck(lesModel.coeffDict().lookupOrDefault<scalar>("Ck",0.094));
    //Info << "Ck = " << Ck << nl;
    Ceta = Ceta0*1e-3/(Ck+1e-20);
    Info << "    Ceta = " << Ceta << nl;
   
    return;
}

Foam::scalar //tmp<Foam::volScalarField>
Foam::radiation::spSootRadFracEmission::interpolatePsiR(const label i, const label j, const label k, const scalar sRate, const scalar pLoss) const
{
	label i1 = i;
	label i2 = i+1;
    	label j1 = j;
    	label j2 = j+1;
	scalar f_x_y1 = ((sr_[i2]-sRate)/(sr_[i2]-sr_[i1]))*psiR(i1,j1,k) + 
         		((sRate-sr_[i1])/(sr_[i2]-sr_[i1]))*psiR(i2,j1,k);
	scalar f_x_y2 = ((sr_[i2]-sRate)/(sr_[i2]-sr_[i1]))*psiR(i1,j2,k) + 
         		((sRate-sr_[i1])/(sr_[i2]-sr_[i1]))*psiR(i2,j2,k); 
	return ((hp_[j2]-pLoss)/(hp_[j2]-hp_[j1]))*f_x_y1 +
        	((pLoss-hp_[j1])/(hp_[j2]-hp_[j1]))*f_x_y2;
}

//void Foam::radiation::spSootRadFracEmission::computeQr() const
Foam::tmp<Foam::volScalarField>
Foam::radiation::spSootRadFracEmission::computeQr() const
{
    tmp<volScalarField> PsiR
    (
        new volScalarField
        (
            IOobject
            (
                "PsiR",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    //volScalarField PsiR
    //(
    //    IOobject
    //    (
    //        "PsiR",
    //        mesh().time().timeName(),
    //        mesh(),
    //        IOobject::NO_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    mesh(),
    //    dimensionedScalar("zero", dimless, 0.0)
    //);
    volScalarField priorLoss
    (
        IOobject
        (
            "priorLoss",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    );
    // volScalarField strainRate
    // (
    //     IOobject
    //     (
    //         "strainRate",
    //         mesh().time().timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedScalar("zero", pow(dimTime,-1), 1e3)
    // );

    const psiReactionThermo& thermo = mesh_.lookupObject<psiReactionThermo>("thermophysicalProperties");
    //access species mass fractions
    const PtrList<volScalarField>& Y = thermo.composition().Y();
    //access specie thermo data
    const PtrList<gasHThermoPhysics> & specieThermo =
    	dynamic_cast<const reactingMixture<gasHThermoPhysics>&>  (thermo).speciesData();

    //tmp<volScalarField> Thc_(Thermo.hc());
    //const volScalarField& Hc_ = Thc_();

    // get index of O2 in mixture
    label indexO2 = dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["O2"];
    // get index of N2 in mixture
    label indexN2 = dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["N2"];
    // Computing O2 chemical enthalpy
    const volScalarField& O2 = Y[indexO2];
    scalar W_O2 = specieThermo[indexO2].W();
    scalar h_O2 = specieThermo[indexO2].hc();
    scalar hcO2 = h_O2/W_O2;
    // Computing N2 chemical enthalpy
    const volScalarField& N2 = Y[indexN2];
    scalar W_N2 = specieThermo[indexN2].W();
    scalar h_N2 = specieThermo[indexN2].hc();
    scalar hcN2 = h_N2/W_N2;
   
    // Computing mixture fraction
    scalar s = dynamic_cast<const singleStepReactingMixture<gasHThermoPhysics>&> (thermo).s().value();
    // get index of fuel in mixture
    label fuelIndex = dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()[fuel_];
    const volScalarField& Fu = Y[fuelIndex];
    scalar Wu_ = specieThermo[fuelIndex].W();
    scalar Hu_ = specieThermo[fuelIndex].hc();
    scalar hcFuel = Hu_/Wu_;
    volScalarField Ft 
    (
    IOobject
    (
        "Ft",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (Fu*s-O2+YO2Inf)/(s+YO2Inf)
    );
    Ft.max(0.0);
    Ft.min(1.0);

    // Computing strain rate
    const psiThermo& Thermo = mesh_.lookupObject<psiThermo>("thermophysicalProperties");
    tmp<volScalarField> tmu(Thermo.mu());
    const volScalarField& mu_ = tmu();
    tmp<volScalarField> tepsilon(turbulence_.epsilon());
    const volScalarField& epsilon = tepsilon();
    tmp<volScalarField> trho(thermo_.rho());
    const volScalarField& rho_ = trho();
    dimensionedScalar TINY
    (
        "TINY",
        dimMass/(dimLength*dimTime),
        1.0e-9
    );
    strainRate_ = Ceta*pow(epsilon*rho_/(mu_+TINY),0.5);

    // Computing prior loss fraction
    tmp<volScalarField> thc(thermo_.hc());
    const volScalarField& hc = thc();
    tmp<volScalarField> the(thermo_.he());
    const volScalarField& he = the();
    dimensionedScalar enthUnits
    (
        "enthUnits",
        dimVelocity*dimVelocity,
        1.0
    );
    dimensionedScalar TINY2
    (
        "TINY2",
        dimVelocity*dimVelocity,
        1.0e1 //-6
    );
    volScalarField correctHc = (O2*hcO2 + N2*hcN2)*enthUnits;
    priorLoss = (he + hc - Ft*hcFuel*enthUnits - correctHc)/(hc - Ft*hcFuel*enthUnits + TINY2 - correctHc);

    PsiR.ref().ref() = 0.0*priorLoss;
    forAll(strainRate_,cellI)
    {
        //strainRate[cellI] = 15;
        //priorLoss[cellI] = 0.05;
        scalar Ydummy = 0.0;
        if ( globalYO2_ ) {
            Ydummy = YO2Inf;
        }
        else {
            scalar Vsum = 0.0;
            Ydummy = O2[cellI]*mesh().V()[cellI];
            forAll(mesh().cellCells()[cellI],nbri)
            {
                Ydummy += O2[nbri]*mesh().V()[nbri];
                Vsum += mesh().V()[nbri];
            }
            Ydummy /= (Vsum + mesh().V()[cellI]); //= O2[cellI];
        }
        //Info << "Ydummy: " << Ydummy << nl;
        if (Ydummy > tableOxyMassFracs_[0] || Ydummy < tableOxyMassFracs_[tableOxyMassFracs_.size()-1]) {
            PsiR.ref().ref()[cellI] = 0.0;
        }
        else {
		if(Ft[cellI] < 1e-6) {
			priorLoss[cellI] = 1e3;
		}
		label i1,j1;
		if(strainRate_[cellI] > sr_[sr_.size()-1]) {
			strainRate_[cellI] = sr_[sr_.size()-1];
			i1 = sr_.size()-2;
		}
		if(strainRate_[cellI] < sr_[0]) {
			strainRate_[cellI] = sr_[0];
			i1 = 0;
		}
		for(label i=0;i<sr_.size();i++)
		{
			i1 = i;
			if(strainRate_[cellI] >= sr_[i1] && strainRate_[cellI] <= sr_[i1+1]) {
				break;
			}
		}
		if(priorLoss[cellI] > hp_[hp_.size()-1]) {
			priorLoss[cellI] = hp_[hp_.size()-1];
			j1 = hp_.size()-2;
		}
		if(priorLoss[cellI] < hp_[0]) {
			priorLoss[cellI] = hp_[0];
			j1 = 0;
		}
		for(label j=0;j<hp_.size();j++)
		{
			j1 = j;
			if(priorLoss[cellI] >= hp_[j1] && priorLoss[cellI] <= hp_[j1+1]) {
				break;
			}
		}

		forAll(tableOxyMassFracs_,k)
		{
			if (Ydummy >= tableOxyMassFracs_[k+1] && Ydummy < tableOxyMassFracs_[k]) {
				scalar weight = (Ydummy - tableOxyMassFracs_[k+1])/(tableOxyMassFracs_[k] - tableOxyMassFracs_[k+1]);
				PsiR.ref().ref()[cellI] = weight*interpolatePsiR(i1,j1,k,strainRate_[cellI],priorLoss[cellI]) +
					(1-weight)*interpolatePsiR(i1,j1,k+1,strainRate_[cellI],priorLoss[cellI]);
				break;
			}
			else if (Ydummy == tableOxyMassFracs_[k]) {
				PsiR.ref().ref()[cellI] = interpolatePsiR(i1,j1,k,strainRate_[cellI],priorLoss[cellI]);
				break;
			}
		}
	}
	if(PsiR.ref().ref()[cellI] < 0.0 || PsiR.ref().ref()[cellI] > 1.0) {
		PsiR.ref().ref()[cellI] = 0.0;
	}
	if(priorLoss[cellI] >= 1.0 || priorLoss[cellI] <= -1){
		priorLoss[cellI] = 0.0;
	}
    }

    // Debug output
    static bool DEBUG = false;
    if ( DEBUG ) {
    	Info << "Y_O2: " << O2 << nl;
    	Info << "Y_N2: " << N2 << nl;
    	//Info << "Y: " << Y << nl;
    	Info << "Y_F: " << Fu << nl;
    	Info << "s: " << s << nl;
    	Info << "Wu_: " << Wu_ << nl;
    	Info << "Hu_: " << Hu_ << nl;
    	Info << "YO2Inf: " << YO2Inf << nl;
    	Info << "epsilon: " << epsilon << nl;
    	Info << "mu: " << mu_ << nl;
    	Info << "rho: " << rho_ << nl;
    	Info << "Ceta: " << Ceta << nl;
    	Info << "strainRate_: " << strainRate_ << nl;
    	Info << "he: " << he << nl;
    	//Info << "Hc_: " << Hc_ << nl;
    	Info << "hc: " << hc << nl;
    	Info << "Ft: " << Ft << nl;
    	Info << "hcFuel: " << hcFuel << nl;
    	Info << "correctHc: " << correctHc << nl;
    	Info << "priorLoss: " << priorLoss << nl;
    }

if(mesh().time().outputTime()){
	strainRate_.write();
	priorLoss.write();
}

    //const volScalarField& PsiR1 = PsiR;
    //Info << "PsiR: " << PsiR1 << nl;
    return PsiR;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::spSootRadFracEmission::ECont(const label bandI) const
{

//static bool callEcont = true;

    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    tmp<volScalarField> psiR(computeQr());
    const volScalarField& PsiR = psiR();
    if(mesh().time().outputTime()){
	PsiR.write();
    }
    scalar RadFraction = 0;

    if (radScaling)
    {

        //TODO: this doesn't need to be recomputed for each ILambda solve

        const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");

        scalar mlr1(0.0);

        // kvm added this ....

        forAll(patchName1_,i)
        {
            const label patchI = mesh_.boundaryMesh().findPatchID(patchName1_[i]);
            if(patchI<0)
            {
                FatalErrorIn("radScaling.H")
                    << "patch " << patchName1_[i] << " not found" << nl
                    << abort(FatalError);
            }
            mlr1 += -gSum(phi.boundaryField()[patchI]);
        }

        scalar mlr2(0.0);

        forAll(patchName2_,i)
        {
            const label patchI = mesh_.boundaryMesh().findPatchID(patchName2_[i]);
            if(patchI<0)
            {
                FatalErrorIn("radScaling.H")
                    << "patch " << patchName2_[i] << " not found" << nl
                    << abort(FatalError);
            }
            mlr2 += -gSum(phi.boundaryField()[patchI]);
        }

        if(debug)
        {
            Info << "mlr for patches " << patchName1_ << " is " << mlr1 << endl;
            Info << "mlr for patches " << patchName2_ << " is " << mlr2 << endl;
        }

        scalar minRadFrac = min(Ehrr1_,Ehrr2_);

        RadFraction = (mlr1*Ehrr1_ + mlr2*Ehrr2_)
                    / max(SMALL, (mlr1 + mlr2));
        RadFraction = max(minRadFrac,RadFraction);
        //debug Info << "RadFraction " << RadFraction << endl;
    }
    else
    {
        RadFraction = EhrrCoeff_;
    }


    if (mesh_.foundObject<volScalarField>("Qdot"))
    {
        const volScalarField& Qdot =
            mesh_.lookupObject<volScalarField>("Qdot");
        if (Qdot.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            E.ref().ref() = PsiR*Qdot;
            //E.ref().ref() = RadFraction*Qdot;
        }
        else
        {
            Info << "Qdot dimensions incorrect" << endl;
        }

        static word timeName = "null";
        if (timeName != mesh().time().timeName())
        {
    const volScalarField& e = E(); 
    const scalarField& v = mesh().V();
    //volScalarField rad(e*v());
	    
    dimensionedScalar TINY
    (
        "TINY",
        dimEnergy/dimTime/dimVolume,
        1.0e-9
    );
    volScalarField qdot = Qdot+TINY; 
	    
            Info << "Radiant Fraction is " << gSum(e*v)/gSum(qdot*v) << nl;
            //Info << "Radiant Fraction is " << RadFraction << endl;
            timeName = mesh().time().timeName();
        }
    }
    return E;
}


// ************************************************************************* //

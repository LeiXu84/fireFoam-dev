/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "reactingOneDimRPvol.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvMatrices.H"
#include "absorptionEmissionModel.H"
#include "fvcLaplacian.H"
#include "surfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactingOneDimRPvol, 0);

addToRunTimeSelectionTable(pyrolysisModel, reactingOneDimRPvol, mesh);
addToRunTimeSelectionTable(pyrolysisModel, reactingOneDimRPvol, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void reactingOneDimRPvol::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;
    time().controlDict().lookup("maxDi") >> maxDiff_;

    coeffs().lookup("minimumDelta") >> minimumDelta_;

    coeffs().lookup("gasHSource") >> gasHSource_;
    coeffs().lookup("qrHSource") >> qrHSource_;
    useChemistrySolvers_ =
        coeffs().lookupOrDefault<bool>("useChemistrySolvers", true);

    coeffs().lookup("Tcrt") >> Tcrt_;
    coeffs().lookup("TcrtLocal") >> TcrtLocal_;
    coeffs().lookup("deltTLocal") >> deltTLocal_;
    paperToFuelRatio_ = coeffs().lookupOrDefault<scalar>("paperToFuelRatio",0.9);
}


bool reactingOneDimRPvol::read()
{
    if (pyrolysisModel::read())
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


bool reactingOneDimRPvol::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


void reactingOneDimRPvol::updateqr()
{
    // Update local qr from coupled qr field
    qr_ == dimensionedScalar("zero", qr_.dimensions(), 0.0);

    // Retrieve field from coupled region using mapped boundary conditions
    qr_.correctBoundaryConditions();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        scalarField& qrp = qr_.boundaryFieldRef()[patchi];

        // qr is positive going in the solid
        // If the surface is emitting the radiative flux is set to zero
        qrp = max(qrp, scalar(0.0));
    }

    const vectorField& cellC = regionMesh().cellCentres();

    tmp<volScalarField> kappa = kappaRad();

    // Propagate qr through 1-D regions
    label localPyrolysisFacei = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        const scalarField& qrp = qr_.boundaryField()[patchi];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchi];

        forAll(qrp, facei)
        {
            const scalar qr0 = qrp[facei];
            point Cf0 = Cf[facei];
            const labelList& cells = boundaryFaceCells_[localPyrolysisFacei++];
            scalar kappaInt = 0.0;
            forAll(cells, k)
            {
                const label celli = cells[k];
                const point& Cf1 = cellC[celli];
                const scalar delta = mag(Cf1 - Cf0);
                kappaInt += kappa()[celli]*delta;
                qr_[celli] = qr0*exp(-kappaInt);
                Cf0 = Cf1;
            }
        }
    }
}


void reactingOneDimRPvol::updatePhiGas()
{
    phiHsGas_ ==  dimensionedScalar("zero", phiHsGas_.dimensions(), 0.0);
    phiGas_ == dimensionedScalar("zero", phiGas_.dimensions(), 0.0);

    const speciesTable& gasTable = solidChemistry_->gasTable();

    scalar mlossAll = 0;

    forAll(gasTable, gasI)
    {
        tmp<volScalarField> tHsiGas =
            solidChemistry_->gasHs(solidThermo_.p(), solidThermo_.T(), gasI);

        const volScalarField& HsiGas = tHsiGas();

        const volScalarField::Internal& RRiGas =
            solidChemistry_->RRg(gasI);

        label totalFaceId = 0;
        forAll(intCoupledPatchIDs_, i)
        {
            const label patchi = intCoupledPatchIDs_[i];

            scalarField& phiGasp = phiGas_.boundaryFieldRef()[patchi];
            const scalarField& cellVol = regionMesh().V();

            forAll(phiGasp, facei)
            {
                const labelList& cells = boundaryFaceCells_[totalFaceId++];
                scalar massInt = 0.0;
                forAllReverse(cells, k)
                {
                    const label celli = cells[k];
                    massInt += RRiGas[celli]*cellVol[celli];
                    phiHsGas_[celli] += massInt*HsiGas[celli];
                }

                phiGasp[facei] += massInt;
                const label cell0 = cells[0];
                const label cell1 = cells[1];
                const label cell2 = cells[2];
                //phiGasp[facei] += 0.0005*Ys_[0][cell0]*virginD_[cell0]/time().deltaTValue();
                //virginD_[cell0] = 0.9995*virginD_[cell0];
                scalar dtime = time().deltaTValue();

                //-Try double smooth method
                // virginD_[cell0] = accumulative fuel to the buffer
                // virginD_[cell1] = smooth virginD_[cell0]
                // virginD_[cell2] = smooth virginD_[cell1]
                virginD_[cell1] = (virginD_[cell1] + virginD_[cell0]*dtime)/(1.0+dtime);
                phiGasp[facei] += (virginD_[cell1] - virginD_[cell2])/(1.0+dtime);
                mlr_[cell0] = (virginD_[cell1] - virginD_[cell2])/(1.0+dtime);
                virginD_[cell2] = (virginD_[cell2] + virginD_[cell1]*dtime)/(1.0+dtime);
                mlossAll += mlr_[cell0];

                //-Try uniform fuel release rate for each page
                // virginD_[cell0] = accumulative fuel to the buffer
                // virginD_[cell1] = unreleased fuel in the buffer
                // virginD_[cell2] = current fuel release rate in mass/time
                //if((virginD_[cell1]-virginD_[cell2]*dtime) > 0.0)
                //{
                //    phiGasp[facei] = virginD_[cell2];
                //    virginD_[cell1] -= virginD_[cell2]*dtime;
                //}else
                //{
                //    phiGasp[facei] = 0;
                //}

                if (debug)
                {
                    Info<< " Gas : " << gasTable[gasI]
                        << " on patch : " << patchi
                        << " mass produced at face(local) : "
                        <<  facei
                        << " is : " << massInt
                        << " [kg/s] " << endl;
                }
            }
        }
        tHsiGas.ref().clear();
    }
    Info<<"dbg-MLR: "<<tab<<time().value()<<tab<<mlossAll<<endl;
}


void reactingOneDimRPvol::updateFields()
{
    if (qrHSource_)
    {
        updateqr();
    }

    updatePhiGas();
}


void reactingOneDimRPvol::updateMesh(const scalarField& mass0)
{
    if (!moveMesh_)
    {
        return;
    }

    const scalarField newV(mass0/rho_);

    Info<< "Initial/final volumes = " << gSum(regionMesh().V()) << ", "
        << gSum(newV) << " [m3]" << endl;

    // move the mesh
    const labelList moveMap = moveMesh(regionMesh().V() - newV, minimumDelta_);

    // flag any cells that have not moved as non-reacting
    forAll(moveMap, i)
    {
        if (moveMap[i] == 0)
        {
            solidChemistry_->setCellReacting(i, false);
        }
    }
}


void reactingOneDimRPvol::solveContinuity()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    const scalarField mass0 = rho_*regionMesh().V();

    fvScalarMatrix rhoEqn
    (
          fvm::ddt(rho_)
        + fvc::div(phiPyrolysis_)
          ==
        - solidChemistry_->RRg()
    );

    if (regionMesh().moving())
    {
        surfaceScalarField phiRhoMesh
        (
            fvc::interpolate(rho_)*regionMesh().phi()
        );

        rhoEqn += fvc::div(phiRhoMesh);
    }

    rhoEqn.solve();

    updateMesh(mass0);
}


void reactingOneDimRPvol::solveSpeciesMass()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho_, Yi)
          + fvc::div(phiPyrolysis_, Yi)
         ==
            solidChemistry_->RRs(i)
        );

        if (regionMesh().moving())
        {
            surfaceScalarField phiYiRhoMesh
            (
                fvc::interpolate(Yi*rho_)*regionMesh().phi()
            );

            YiEqn += fvc::div(phiYiRhoMesh);

        }

        YiEqn.solve(regionMesh().solver("Yi"));
        Yi.max(0.0);
        Yt += Yi;
    }

    Ys_[Ys_.size() - 1] = 1.0 - Yt;

    //-Correct density from mixture density
    forAll(rho_, celli)
    {
        scalar rc = 0;
        forAll(Ys_, i)
        {
            rc += Ys_[i][celli]/solidThermo_.composition().rho(i, 1.0e5, 300.0);
        }
        //Info<<"dbg-rho: "<<rho_[celli]<<tab<<1.0/rc<<endl;
        rho_[celli] = 1.0/rc;
    }
}


void reactingOneDimRPvol::solveEnergy()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    tmp<volScalarField> alpha(solidThermo_.alpha());

    dimensionedScalar Cp0("Cp0", dimEnergy/dimMass/dimTemperature, solidThermo_.composition().Cp(0, 1.0e05, 300.) );
    dimensionedScalar Cp1("Cp1", dimEnergy/dimMass/dimTemperature, solidThermo_.composition().Cp(1, 1.0e05, 300.) );

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho_, h_)
      + fvc::div(phiPyrolysis_, h_)
      - fvm::laplacian(alpha, h_)
      + fvc::laplacian(alpha, h_)
      - fvc::laplacian(kappa(), T())
     ==
        chemistryQdot_
//      - fvm::Sp(solidChemistry_->RRg(), h_)
      + solidChemistry_->RRs(0)*T()*Cp0
      + solidChemistry_->RRs(1)*T()*Cp1
    );

    if (gasHSource_)
    {
        const surfaceScalarField phiGas(fvc::interpolate(phiHsGas_));
        hEqn += fvc::div(phiGas);
    }

    if (qrHSource_)
    {
        const surfaceScalarField phiqr(fvc::interpolate(qr_)*nMagSf());
        hEqn += fvc::div(phiqr);
    }

    if (regionMesh().moving())
    {
        surfaceScalarField phihMesh
        (
            fvc::interpolate(rho_*h_)*regionMesh().phi()
        );

        hEqn += fvc::div(phihMesh);
    }

    hEqn.relax();
    hEqn.solve();
}


void reactingOneDimRPvol::calculateMassTransfer()
{
    totalGasMassFlux_ = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        totalGasMassFlux_ += gSum(phiGas_.boundaryField()[patchi]);
    }

    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistryQdot_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingOneDimRPvol::reactingOneDimRPvol
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, regionType),
    solidChemistry_(basicSolidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    radiation_(radiation::radiationModel::New(solidThermo_.T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_.rho(),
        zeroGradientFvPatchScalarField::typeName
    ),
    Ys_(solidThermo_.composition().Y()),
    h_(solidThermo_.he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    qr_
    (
        IOobject
        (
            "qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
        //dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        //zeroGradientFvPatchVectorField::typeName
    ),
    
    //-Roll paper
    Tcrt_(600),
    TcrtLocal_(500),
    deltTLocal_(5),
    paperToFuelRatio_(0.9),

    Upyrolysis_
    (
        IOobject
        (
            "Upyrolysis",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector("zero", dimLength/dimTime, vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),

    phiPyrolysis_
    (
        IOobject
        (
            "phiPyrolysis",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),
    
    page_
    (
        IOobject
        (
            "page",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimMass, 0.0)
    ),

    blockFactor_
    (
        IOobject
        (
            "blockFactor",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    wd_
    (
        IOobject
        (
            "weightLocal",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Ts0_
    (
        IOobject
        (
            "Ts0",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
    dimensionedScalar("zero", dimTemperature, 1000.0)
    ),
    
    TdLocal_
    (
        IOobject
        (
            "TdLocal",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 1000.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    TdThin_
    (
        IOobject
        (
            "TdThin",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    TdFirst_
    (
        IOobject
        (
            "TdFirst",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    QnetSmooth_
    (
        IOobject
        (
            "QnetSmooth",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    dMoved_
    (
        IOobject
        (
            "dMoved",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimLength, -1.0)
    ),

    virginD_
    (
        IOobject
        (
            "virginD",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mlr_
    (
        IOobject
        (
            "mlr",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
    dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0)
    ),

    Qnet_
    (
        IOobject
        (
            "Qnet",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    gasHSource_(false),
    qrHSource_(false),
    useChemistrySolvers_(true)
{
    if (active_)
    {
        read();
    }
}

reactingOneDimRPvol::reactingOneDimRPvol
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, dict, regionType),
    solidChemistry_(basicSolidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    radiation_(radiation::radiationModel::New(solidThermo_.T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_.rho(),
        zeroGradientFvPatchScalarField::typeName
    ),
    Ys_(solidThermo_.composition().Y()),
    h_(solidThermo_.he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    qr_
    (
        IOobject
        (
            "qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    //-Roll paper
    Tcrt_(600),
    TcrtLocal_(500),
    deltTLocal_(5),
    paperToFuelRatio_(0.9),

    Upyrolysis_
    (
        IOobject
        (
            "Upyrolysis",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedVector("zero", dimLength/dimTime, vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),

    phiPyrolysis_
    (
        IOobject
        (
            "phiPyrolysis",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),
    
    page_
    (
        IOobject
        (
            "page",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    blockFactor_
    (
        IOobject
        (
            "blockFactor",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    wd_
    (
        IOobject
        (
            "weightLocal",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Ts0_
    (
        IOobject
        (
            "Ts0",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
    dimensionedScalar("zero", dimTemperature, 1000.0)
    ),

    TdLocal_
    (
        IOobject
        (
            "TdLocal",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 1000.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    TdThin_
    (
        IOobject
        (
            "TdThin",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    TdFirst_
    (
        IOobject
        (
            "TdFirst",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    QnetSmooth_
    (
        IOobject
        (
            "QnetSmooth",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    dMoved_
    (
        IOobject
        (
            "dMoved",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimLength, -1.0)
    ),

    virginD_
    (
        IOobject
        (
            "virginD",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mlr_
    (
        IOobject
        (
            "mlr",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
    dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0)
    ),

    Qnet_
    (
        IOobject
        (
            "Qnet",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    gasHSource_(false),
    qrHSource_(false),
    useChemistrySolvers_(true)
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactingOneDimRPvol::~reactingOneDimRPvol()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar reactingOneDimRPvol::addMassSources(const label patchi, const label facei)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchi)
        {
            index = i;
            break;
        }
    }

    const label localPatchId =  intCoupledPatchIDs_[index];

    const scalar massAdded = phiGas_.boundaryField()[localPatchId][facei];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


scalar reactingOneDimRPvol::solidRegionDiffNo() const
{
    scalar DiNum = -GREAT;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            sqr(regionMesh().surfaceInterpolation::deltaCoeffs())
           *fvc::interpolate(kappa())
           /fvc::interpolate(Cp()*rho_)
        );

        DiNum = max(KrhoCpbyDelta.primitiveField())*time().deltaTValue();
    }

    return DiNum;
}


scalar reactingOneDimRPvol::maxDiff() const
{
    return maxDiff_;
}


const volScalarField& reactingOneDimRPvol::rho() const
{
    return rho_;
}


const volScalarField& reactingOneDimRPvol::T() const
{
    return solidThermo_.T();
}


const tmp<volScalarField> reactingOneDimRPvol::Cp() const
{
    return solidThermo_.Cp();
}


tmp<volScalarField> reactingOneDimRPvol::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> reactingOneDimRPvol::kappa() const
{
    return solidThermo_.kappa();
}


const surfaceScalarField& reactingOneDimRPvol::phiGas() const
{
    return phiGas_;
}


void reactingOneDimRPvol::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();

    // Initialise all cells as able to react
    forAll(h_, celli)
    {
        solidChemistry_->setCellReacting(celli, true);
    }

    //// Correct Boundary conditions
    //page_.correctBoundaryConditions();
    //blockFactor_.correctBoundaryConditions();
    //wd_.correctBoundaryConditions();
    //QnetSmooth_.correctBoundaryConditions();

    // Move material up
    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();
    label localPyrolysisFaceI = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch pp = bm[patchi];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchi];
        const scalarField& pageNum = page_.boundaryField()[patchi];
        const scalarField& QnetFace = Qnet_.boundaryField()[patchi];

        // Get local weighting from film model
        typedef regionModels::surfaceFilmModels::surfaceFilmModel 
            surfaceFilmModelType;
        const regionModels::regionModel& surFilmRegion =
            db().time().lookupObject<regionModels::regionModel>
            ("surfaceFilmProperties");
        const surfaceFilmModelType& surFilmModel = 
            dynamic_cast<const surfaceFilmModelType&>(surFilmRegion);
        const surfaceFilmModelType& surFilm = 
            const_cast<surfaceFilmModelType&>(surFilmModel);
        scalarList wLocal(pp.faceCells().size(), 0.0);
        wLocal = mapRegionPatchInternalField<scalar>
            (
                surFilm,
                "paperWeight",
                patchi,
                true
            );

    forAll(pp, facei)
    {
        const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI];
        const labelList& faces = boundaryFaceFaces_[localPyrolysisFaceI];
        localPyrolysisFaceI++;
        const vector sf = pp.faceAreas()[facei];
        const label cell0 = cells[0];
        const label cell1 = cells[1];
        const label cell2 = cells[2];

        scalar paperThick = regionMesh().V()[cell0]/mag(sf);
        scalar uMoveMax = paperThick/time().deltaTValue();
        scalar uMove(0.201*uMoveMax);   //Finish moving in 5 time steps
        scalar Tsurface = T()[cell0];
        scalar dTdtS = (Tsurface - Ts0_[cell0])/time().deltaTValue();
        Ts0_[cell0] = Tsurface;
        Qnet_[cell0] = QnetFace[facei];
        QnetSmooth_[cell0] = (QnetSmooth_[cell0] + Qnet_[cell0]*time().deltaTValue())/(1.0+time().deltaTValue());
        //scalar fa0 = 1.0 - 30.0/(39.0+page_[cell0]);
        //if(page_[cell0] == 0)
        //{
        //    fa0 = 0;
        //}
        //scalar tshrink = 0.5;

        scalar fa0 = 0.5 - 10.0/(20.0+page_[cell0]);    // Ning
        //scalar fa0 = 1.0 - 35.0/(39.0+page_[cell0]);  // Dong
        //scalar fa0 = 0;                   // No Page Blocking

        scalar tshrink = paperThick*1.5e8/(1.0e4+fabs(QnetSmooth_[cell0]));
        if(tshrink>5.0){Info<<"Tau too large! "<<tshrink<<endl; tshrink = 5.0;}
        if(tshrink<0.05){Info<<"Tau too small! "<<tshrink<<endl; tshrink = 0.05;}

        bool updateLocal = false;

        //-Local delamination temperature
        //scalar zz = Cf[facei].z();
        //if(zz < 0.0) Info<<"Error! "<<zz<<endl;
        ////wLocal[facei] = wLocal[facei]*sqrt(20.0*zz)/5.0;
        //-Discontinuity between tiers
        //if( ((zz>2.1)&&(zz<2.13)) || ((zz>4.2)&&(zz<4.23)) )
        //{
        //    wLocal[facei] = 0;
        //}
        if(wd_[cell0] != wLocal[facei])
        {
        updateLocal = true;
        wd_[cell0] = wLocal[facei];
        TdLocal_[cell0] = Tsurface + (Tcrt_-Tsurface)/(1.0+wd_[cell0]);
        Info<<"dbg-local "<<time().value()<<tab
            <<Tsurface<<tab<<TdLocal_[cell0]<<tab<<dTdtS<<tab
            <<wd_[cell0]<<tab<<facei<<endl;
        }

        scalar thermoInertia = rho_[cell0]*paperThick*2500;
        if(wd_[cell0] > 0)
        {
        TdThin_[cell0] += QnetSmooth_[cell0]*time().deltaTValue()/thermoInertia;
        }
        else
        {
        TdThin_[cell0] = T()[cell0];
        }

        if(dMoved_[cell0] < 0)
        {
        bool delaminate = false;
        if(Tsurface > Tcrt_ || Tsurface > TdLocal_[cell0])  // 2-D Model
        //if(Tsurface > Tcrt_)                  // 1-D Model
        {
            delaminate = true;
            if(page_[cell0] == 0){TdFirst_[cell0] = Tsurface;}
        }
        if(delaminate)
        {
            forAll(cells, cI)
            {
            Upyrolysis_[cells[cI]] = uMove*sf/mag(sf);
            }
            //Info<<"dbg-dela: "<<time().value()<<tab<<localPyrolysisFaceI<<tab
            //<<page_[cell0]<<tab<<blockFactor_[cell0]<<tab<<fa0<<tab<<tshrink<<endl;
            dMoved_[cell0] = 0;
            scalar pM=rho_[cell0]*regionMesh().V()[cell0]*Ys_[0][cell0]*paperToFuelRatio_;
            virginD_[cell0] += pM;
            page_[cell0]++;
            blockFactor_[cell0] = 1.0;
            //Info<<"dbg-da: "<<time().value()<<tab<<localPyrolysisFaceI<<tab
            //<<page_[cell0]<<tab<<blockFactor_[cell0]<<tab<<fa0<<tab<<tshrink<<endl;
        }
        }else if(dMoved_[cell0] > paperThick)
        {
            forAll(cells, cI)
            {
           Upyrolysis_[cells[cI]] = vector::zero;
            }
            dMoved_[cell0] = -1;
        }else
        {
        dMoved_[cell0] += uMove*time().deltaTValue();
        }
        blockFactor_[cell0] += - (blockFactor_[cell0]-fa0)*time().deltaTValue()/tshrink;
        //blockFactor_[cell0] = max(0.0, fa0);

        if(page_[cell0] == 0)
        {
            blockFactor_[cell0] = 0;
        }
    }   //for all faces on bm.
    }   //for all patches of pyrolysis zone.

    // Correct Boundary conditions
    page_.correctBoundaryConditions();
    blockFactor_.correctBoundaryConditions();
    wd_.correctBoundaryConditions();
    TdThin_.correctBoundaryConditions();
    TdFirst_.correctBoundaryConditions();
    QnetSmooth_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    Upyrolysis_.correctBoundaryConditions();
    phiPyrolysis_ = linearInterpolate(rho_*Upyrolysis_) & regionMesh().Sf();
    //Info<<"PhiPyrolysis: \n"<<phiPyrolysis_<<endl;
}


void reactingOneDimRPvol::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    if (useChemistrySolvers_)
    {
        solidChemistry_->solve(time().deltaTValue());
    }
    else
    {
        solidChemistry_->calculate();
    }

    solveContinuity();

    chemistryQdot_ = solidChemistry_->Qdot()();

    updateFields();

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }

    calculateMassTransfer();

    solidThermo_.correct();

    Info<< "pyrolysis min/max(T) = "
        << min(solidThermo_.T().internalField())
        << ", "
        << max(solidThermo_.T().internalField())
        << endl;
}


void reactingOneDimRPvol::info()
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << endl;

    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //

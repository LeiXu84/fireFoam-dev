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

#include "reactingOneDimRollPaper.H"
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

#define DEBUG(x) {                                              \
            std::streamsize p = std::cout.precision();              \
            std::ios::fmtflags myFlags;                             \
            myFlags = cout.flags();                                 \
            std::cout.precision(10);                                \
            std::cout.setf(std::ios::fixed,std::ios::floatfield);   \
            std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
            std::cout << "p" << Pstream::myProcNo();                \
            std::cout << " " << #x " = " << x << std::endl;         \
            std::cout.precision(p);                                 \
            std::cout.flags(myFlags);                               \
        }
#define TRACE(s) {                                              \
            std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
            std::cout << "p" << Pstream::myProcNo();                \
            std::cout << " " << #s << std::endl;                    \
            s;                                                      \
        }

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactingOneDimRollPaper, 0);

addToRunTimeSelectionTable(pyrolysisModel, reactingOneDimRollPaper, mesh);
addToRunTimeSelectionTable(pyrolysisModel, reactingOneDimRollPaper, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void reactingOneDimRollPaper::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;
    time().controlDict().lookup("maxDi") >> maxDiff_;

    coeffs().lookup("minimumDelta") >> minimumDelta_;

    coeffs().lookup("gasHSource") >> gasHSource_;
    coeffs().lookup("QrHSource") >> QrHSource_;
    useChemistrySolvers_ =
        coeffs().lookupOrDefault<bool>("useChemistrySolvers", true);

    coeffs().lookup("Tcrt") >> TcriticalDelamination_;
    paperToFuelRatio_ = coeffs().lookupOrDefault<scalar>("paperToFuelRatio",0.9);
    Hpyrolysis_ = coeffs().lookupOrDefault<scalar>("Hpyrolysis",1e6);
    blockOfBurningPaper_ = coeffs().lookupOrDefault<scalar>("blocking",0.5);
}


bool reactingOneDimRollPaper::read()
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


bool reactingOneDimRollPaper::read(const dictionary& dict)
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


void reactingOneDimRollPaper::updateQr()
{
    // Update local Qr from coupled Qr field
    Qr_ == dimensionedScalar("zero", Qr_.dimensions(), 0.0);

    // Retrieve field from coupled region using mapped boundary conditions
    Qr_.correctBoundaryConditions();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];

        scalarField& Qrp = Qr_.boundaryFieldRef()[patchI];

        // Qr is positive going in the solid
        // If the surface is emitting the radiative flux is set to zero
        Qrp = max(Qrp, scalar(0.0));
    }

    const vectorField& cellC = regionMesh().cellCentres();

    tmp<volScalarField> kappa = kappaRad();

    // Propagate Qr through 1-D regions
    label localPyrolysisFaceI = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];

        const scalarField& Qrp = Qr_.boundaryField()[patchI];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchI];

        forAll(Qrp, faceI)
        {
            const scalar Qr0 = Qrp[faceI];
            point Cf0 = Cf[faceI];
            const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI++];
            scalar kappaInt = 0.0;
            forAll(cells, k)
            {
                const label cellI = cells[k];
                const point& Cf1 = cellC[cellI];
                const scalar delta = mag(Cf1 - Cf0);
                kappaInt += kappa()[cellI]*delta;
                Qr_[cellI] = Qr0*exp(-kappaInt);
                Cf0 = Cf1;
            }
        }
    }
}


void reactingOneDimRollPaper::updatePhiGas()
{
    phiHsGas_ ==  dimensionedScalar("zero", phiHsGas_.dimensions(), 0.0);
    phiGas_ == dimensionedScalar("zero", phiGas_.dimensions(), 0.0);

    const speciesTable& gasTable = solidChemistry_->gasTable();

    scalar mlossAll = 0;

    scalar dtime = time().deltaT().value();
    scalar dtime0 = time().deltaT0().value();

    //scalar relax = 1.0;
    //deltat_ = relax*dtime+(1.0-relax)*deltat_; 
    
    // Info<<"reactingOneDimRollPaper::deltat_: " << tab << time().timeName() << tab << deltat_ << tab << dtime <<endl;

    forAll(gasTable, gasI)
    {
        tmp<volScalarField> tHsiGas =
            solidChemistry_->gasHs(solidThermo_.p(), solidThermo_.T(), gasI);

        const volScalarField& HsiGas = tHsiGas();

        const DimensionedField<scalar, volMesh>& RRiGas =
            solidChemistry_->RRg(gasI);

        label totalFaceId = 0;
        scalar massIntMax = 0.0;
        scalar phiGaspMax = 0.0;
        scalar virginDMax = 0.0;
        forAll(intCoupledPatchIDs_, i)
        {
            const label patchI = intCoupledPatchIDs_[i];

            scalarField& phiGasp = phiGas_.boundaryFieldRef()[patchI];
            const scalarField& cellVol = regionMesh().V();

            forAll(phiGasp, faceI)
            {
                const labelList& cells = boundaryFaceCells_[totalFaceId++];
                scalar massInt = 0.0;
                forAllReverse(cells, k)
                {
                    const label cellI = cells[k];
                    massInt += RRiGas[cellI]*cellVol[cellI];
                    phiHsGas_[cellI] += massInt*HsiGas[cellI];
                }
                massIntMax = max(massIntMax,massInt);

                phiGasp[faceI] += massInt;
                const label cell0 = cells[0];
		phiGasp[faceI] = massReleaseRate_[cell0];
		phiGaspMax = max(phiGaspMax, phiGasp[faceI]);

                if (debug)
                {
                    Info<< " Gas : " << gasTable[gasI]
                        << " on patch : " << patchI
                        << " mass produced at face(local) : "
                        <<  faceI
                        << " is : " << massInt
                        << " [kg/s] " << endl;
                }
            }
        }
        tHsiGas.ref().clear();
 
        reduce(massIntMax,maxOp<scalar>());
        reduce(phiGaspMax,maxOp<scalar>());
        Info<<"reactingOneDimRollPaper::massIntMax[" << gasI << "]: " << tab << time().timeName() << tab << massIntMax <<endl;
        Info<<"reactingOneDimRollPaper::phiGaspMax[" << gasI << "]: " << tab << time().timeName() << tab << phiGaspMax <<endl;
        Info<<"reactingOneDimRollPaper::virginDMax[" << gasI << "]: " << tab << time().timeName() << tab << virginDMax <<endl;
    }
}


void reactingOneDimRollPaper::updateFields()
{
    if (QrHSource_)
    {
        updateQr();
    }

    updatePhiGas();
}


void reactingOneDimRollPaper::updateMesh(const scalarField& mass0)
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


void reactingOneDimRollPaper::solveContinuity()
{
    if (debug)
    {
        Info<< "reactingOneDimRollPaper::solveContinuity()" << endl;
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


void reactingOneDimRollPaper::solveSpeciesMass()
{
    if (debug)
    {
        Info<< "reactingOneDimRollPaper::solveSpeciesMass()" << endl;
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
    forAll(rho_, cellI)
    {
        scalar rc = 0;
        forAll(Ys_, i)
        {
            rc += Ys_[i][cellI]/solidThermo_.composition().rho(i, 1.0e5, 300.0);
        }
        //Info<<"dbg-rho: "<<rho_[cellI]<<tab<<1.0/rc<<endl;
        rho_[cellI] = 1.0/rc;
    }
}


void reactingOneDimRollPaper::solveEnergy()
{
    if (debug)
    {
        Info<< "reactingOneDimRollPaper::solveEnergy()" << endl;
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
        chemistrySh_
//      - fvm::Sp(solidChemistry_->RRg(), h_)
      + solidChemistry_->RRs(0)*T()*Cp0
      + solidChemistry_->RRs(1)*T()*Cp1
    );

    if (gasHSource_)
    {
        const surfaceScalarField phiGas(fvc::interpolate(phiHsGas_));
        hEqn += fvc::div(phiGas);
    }

    if (QrHSource_)
    {
        const surfaceScalarField phiQr(fvc::interpolate(Qr_)*nMagSf());
        hEqn += fvc::div(phiQr);
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


void reactingOneDimRollPaper::calculateMassTransfer()
{
    totalGasMassFlux_ = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        totalGasMassFlux_ += gSum(phiGas_.boundaryField()[patchI]);
    }

    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistrySh_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingOneDimRollPaper::reactingOneDimRollPaper
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

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    Qr_
    (
        IOobject
        (
            "Qr",
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
    paperToFuelRatio_(0.9),
    TcriticalDelamination_(600),
    Hpyrolysis_(1e6),
    blockOfBurningPaper_(0.5),

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
    
    paperID_
    (
        IOobject
        (
            "paperID_",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimless, 0.0)
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

    TDelVirginPaper_
    (
        IOobject
        (
            "TDelVirginPaper",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimTemperature, 300.0),
            zeroGradientFvPatchScalarField::typeName
    ),

    Tsurface_
    (
        IOobject
        (
            "Tsurface",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimTemperature, 300.0),
            zeroGradientFvPatchScalarField::typeName
    ),

    movedPaperDistance_
    (
        IOobject
        (
            "movedPaperDistance",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimLength, -1.0)
    ),

    massDelaminatedVirginPaper_
    (
        IOobject
        (
            "massDelaminatedVirginPaper",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimMass, 0.0)
    ),

    massBurningPaper_
    (
        IOobject
        (
            "massBurningPaper",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimMass, 0.0)
    ),

    massReleaseRate_
    (
        IOobject
        (
            "massReleaseRate",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    Qnet_
    (
        IOobject
        (
            "Qnet",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
    ),

    QnetMovingAverage_
    (
        IOobject
        (
            "QnetMovingAverage",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
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
    QrHSource_(false),
    useChemistrySolvers_(true)
{
    if (active_)
    {
        read();
    }
}

reactingOneDimRollPaper::reactingOneDimRollPaper
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

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    Qr_
    (
        IOobject
        (
            "Qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    //-Roll paper
    paperToFuelRatio_(0.9),
    TcriticalDelamination_(600),
    Hpyrolysis_(1e6),
    blockOfBurningPaper_(0.5),

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

    paperID_
    (
        IOobject
        (
            "paperID",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimless, 0.0),
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

    TDelVirginPaper_
    (
        IOobject
        (
            "TDelVirginPaper",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimTemperature, 300.0),
            zeroGradientFvPatchScalarField::typeName
    ),

    Tsurface_
    (
        IOobject
        (
            "Tsurface",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimTemperature, 300.0),
            zeroGradientFvPatchScalarField::typeName
    ),

    movedPaperDistance_
    (
        IOobject
        (
            "movedPaperDistance",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimLength, -1.0)
    ),

    massDelaminatedVirginPaper_
    (
        IOobject
        (
            "massDelaminatedVirginPaper",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimMass, 0.0)
    ),

    massBurningPaper_
    (
        IOobject
        (
            "massBurningPaper",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimMass, 0.0)
    ),

    massReleaseRate_
    (
        IOobject
        (
            "massReleaseRate",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
            dimensionedScalar("zero", dimMass, 0.0)
    ),

    Qnet_
    (
        IOobject
        (
            "Qnet",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
    ),

    QnetMovingAverage_
    (
        IOobject
        (
            "QnetMovingAverage",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
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
    QrHSource_(false),
    useChemistrySolvers_(true)
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactingOneDimRollPaper::~reactingOneDimRollPaper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar reactingOneDimRollPaper::addMassSources(const label patchI, const label faceI)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchI)
        {
            index = i;
            break;
        }
    }

    const label localPatchId =  intCoupledPatchIDs_[index];

    const scalar massAdded = phiGas_.boundaryField()[localPatchId][faceI];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


scalar reactingOneDimRollPaper::solidRegionDiffNo() const
{
    scalar DiNum = -GREAT;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            regionMesh().surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(kappa())
          / fvc::interpolate(Cp()*rho_)
        );

        DiNum = max(KrhoCpbyDelta.primitiveField())*time().deltaTValue();
    }

    return DiNum;
}


scalar reactingOneDimRollPaper::maxDiff() const
{
    return maxDiff_;
}


const volScalarField& reactingOneDimRollPaper::rho() const
{
    return rho_;
}


const volScalarField& reactingOneDimRollPaper::T() const
{
    return solidThermo_.T();
}


const tmp<volScalarField> reactingOneDimRollPaper::Cp() const
{
    return solidThermo_.Cp();
}


tmp<volScalarField> reactingOneDimRollPaper::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> reactingOneDimRollPaper::kappa() const
{
    return solidThermo_.kappa();
}


const surfaceScalarField& reactingOneDimRollPaper::phiGas() const
{
    return phiGas_;
}


void reactingOneDimRollPaper::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();

    // Initialise all cells as able to react
    forAll(h_, cellI)
    {
        solidChemistry_->setCellReacting(cellI, true);
    }

    updateRollPaper();
}


void reactingOneDimRollPaper::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    if (useChemistrySolvers_)
    {
        solidChemistry_->solve
        (
            time().deltaTValue()
        );
    }
    else
    {
        solidChemistry_->calculate();
    }

    solveContinuity();

    chemistrySh_ = solidChemistry_->Sh()();

    updateFields();

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }

    calculateMassTransfer();

    Info<< "dbg-pyrolysisT"<<endl;
    Info<< "pyrolysis min/max(T) before correct = "
        << gMin(solidThermo_.T().internalField())
        << ", "
        << gMax(solidThermo_.T().internalField())
        << endl;

    solidThermo_.correct();

    Info<< "dbg-pyrolysisTcorrect"<<endl;
    Info<< "pyrolysis min/max(T) = "
        << gMin(solidThermo_.T().internalField())
        << ", "
        << gMax(solidThermo_.T().internalField())
        << endl;
}

void reactingOneDimRollPaper::updateRollPaper()
{
    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();
    label localPyrolysisFaceI = 0;

    scalar QnetMax = -VGREAT;
    scalar QnetMin =  VGREAT;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        const polyPatch pp = bm[patchI];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchI];
        const scalarField& QnetFace = Qnet_.boundaryField()[patchI];

        typedef regionModels::surfaceFilmModels::surfaceFilmModel 
	    surfaceFilmModelType;
        const regionModels::regionModel& surFilmRegion =
            db().time().lookupObject<regionModels::regionModel>
            ("surfaceFilmProperties");
        const surfaceFilmModelType& surFilmModel = 
            dynamic_cast<const surfaceFilmModelType&>(surFilmRegion);
        const surfaceFilmModelType& surFilm = 
            const_cast<surfaceFilmModelType&>(surFilmModel);

        scalarList latestThinArea(pp.faceCells().size(), 0.0);    //Thermally thin zone
        latestThinArea = mapRegionPatchInternalField<scalar>
                         (
                             surFilm,
                             "peelingZone",
                             patchI,
                             true
                         );

        scalarList wetArea(pp.faceCells().size(), 0.0);    //Thermally thin zone
        wetArea = mapRegionPatchInternalField<scalar>
                  (
                      surFilm,
                      "deltaf",
                      patchI,
                      true
                  );

        forAll(pp, faceI)
        {
            const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI];
            const labelList& faces = boundaryFaceFaces_[localPyrolysisFaceI];
            localPyrolysisFaceI++;
            const vector sf = pp.faceAreas()[faceI];
            const label cell0 = cells[0];
        
            scalar paperThickness = regionMesh().V()[cell0]/mag(sf);
            scalar cellMass = regionMesh().V()[cell0]*rho_[cell0];
            scalar cellMassCombustible = cellMass*paperToFuelRatio_;
            scalar paperThermoInertia = paperThickness
                *rho_[cell0]*solidThermo_.composition().Cp(0, 1.0e05, 300.);

            scalar uMoveMax = paperThickness/time().deltaTValue();
            scalar uMove(0.201*uMoveMax);    // Finish moving in 5 time steps

	    // Get thermally-thick surface temperature 
            scalar TsurfaceThick = T()[cell0];

            Qnet_[cell0] = QnetFace[faceI];
	    Qnet_[cell0] = min(Qnet_[cell0], 2e5);
            QnetMax = max(QnetMax,Qnet_[cell0]);
            QnetMin = min(QnetMin,Qnet_[cell0]);
	    scalar alphaQ = time().deltaTValue()/(3.0 + time().deltaTValue());
	    QnetMovingAverage_[cell0] = alphaQ*Qnet_[cell0]
	        + (1-alphaQ)*QnetMovingAverage_[cell0];
	    
	    bool readyToMove(false);
	    bool paperMoving(false);
	    bool paperMoved(false);

	    if(movedPaperDistance_[cell0] < 0)
	    {
		readyToMove = true;
		paperMoving = false;
		paperMoved  = false;
	    }
	    else if(movedPaperDistance_[cell0] > paperThickness)
	    {
		readyToMove = false;
		paperMoving = false;
		paperMoved  = true;
	    }
	    else
	    {
		readyToMove = false;
		paperMoving = true;
		paperMoved  = false;
	    }

	    if(paperMoving)
	    {
		movedPaperDistance_[cell0] += uMove*time().deltaTValue();
	    }

	    if(paperMoved)
	    {
		forAll(cells, cI)
		{
		    Upyrolysis_[cells[cI]] = vector::zero;
		}
		movedPaperDistance_[cell0] = -1;
	    }

	    if(readyToMove)
	    {
		bool del_1D(false); // Thermally-driven delamination
		bool del_2D(false); // Paper peeling
		bool haveVirginPaper(false);
		bool haveBurningPaper(false);

		if(massDelaminatedVirginPaper_[cell0] > 0)
		{
		    haveVirginPaper = true;
		}

		if(massBurningPaper_[cell0] > 0)
		{
		    haveBurningPaper = true;
		}

		if(haveVirginPaper) // Heat up vrigin thin paper until ignition
		{
		    //scalar qpaper(qExtra[faceI]+QnetMovingAverage_[cell0]);
		    TDelVirginPaper_[cell0] += (1.0-blockFactor_[cell0])*time().deltaTValue()
		        *Qnet_[cell0]/paperThermoInertia;
		    TDelVirginPaper_[cell0] = max(TDelVirginPaper_[cell0],273.15);
		    if(TDelVirginPaper_[cell0] >= (TcriticalDelamination_ + 1.0))    // Paper is ignited
		    {
			paperID_[cell0] ++;
			massBurningPaper_[cell0] += massDelaminatedVirginPaper_[cell0];
			massDelaminatedVirginPaper_[cell0] = 0;
		    }
		}
		else	// Allow delamination if no virgin thermall-thin paper exists
		{
		    TDelVirginPaper_[cell0] = TsurfaceThick;

		    if(TsurfaceThick > TcriticalDelamination_)
		    {
		        del_1D = true;
		    }

		    if(latestThinArea[faceI] > 0.5)
		    {
		        del_2D = true;
		    }

		    if(del_1D || del_2D)
		    {
			forAll(cells, cI)
			{
			    Upyrolysis_[cells[cI]] = uMove*sf/mag(sf);
			}
			movedPaperDistance_[cell0] = 0;
			massDelaminatedVirginPaper_[cell0] = cellMassCombustible;
		    }
		}

		if(haveBurningPaper)
		{
		    Tsurface_[cell0] = TcriticalDelamination_;
		    //blockFactor_[cell0] = blockOfBurningPaper_;
		    scalar paperLeft = max(0.0, massBurningPaper_[cell0]/cellMassCombustible);
		    blockFactor_[cell0] = min(1.0, paperLeft);
		    massReleaseRate_[cell0] = max(0.0, Qnet_[cell0]*mag(sf)/Hpyrolysis_);
		    scalar dMass = massReleaseRate_[cell0]*time().deltaTValue();
		    if((massBurningPaper_[cell0]-dMass) > 0)
		    {
			massBurningPaper_[cell0] -= dMass;
		    }
		    else    // Paper burns out
		    {
			massBurningPaper_[cell0] = 0;
			massReleaseRate_[cell0] = 0;
		    }
		}
		else
		{
		    blockFactor_[cell0] = 0;
		    if(haveVirginPaper)
		    {
			Tsurface_[cell0] = TDelVirginPaper_[cell0];
		    }
		    else
		    {
			Tsurface_[cell0] = TsurfaceThick;
		    }
		}
	    }
	}
    }    //for all patches of pyrolysis zone.

    reduce(QnetMax,maxOp<scalar>());
    reduce(QnetMin,minOp<scalar>());
    Info << "reactingOneDimRollPaper::QnetMax:" << tab << time().timeName() << tab << QnetMax << endl;
    Info << "reactingOneDimRollPaper::QnetMin:" << tab << time().timeName() << tab << QnetMin << endl;

    // Correct Boundary conditions
    paperID_.correctBoundaryConditions();
    blockFactor_.correctBoundaryConditions();
    Tsurface_.correctBoundaryConditions();
    TDelVirginPaper_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    Upyrolysis_.correctBoundaryConditions();
    phiPyrolysis_ = linearInterpolate(rho_*Upyrolysis_) & regionMesh().Sf();
    //Info<<"PhiPyrolysis: \n"<<phiPyrolysis_<<endl;
}

void reactingOneDimRollPaper::info() const
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

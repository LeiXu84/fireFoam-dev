/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "fvDOM.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOM, 0);
        addToRadiationRunTimeSelectionTables(fvDOM);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOM::initialise()
{
    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        nRay_ = 4*nPhi_*nTheta_;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi/(2.0*nPhi_);
        scalar deltaTheta = pi/nTheta_;
        label i = 0;
        for (label n = 1; n <= nTheta_; n++)
        {
            scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0; // ankur
            // ankur
            if (uniSolidAngles_)
            {
              scalar muBeg = 1. - (n - 1.)*2./nTheta_;   // mu is cos(theta)
              scalar muEnd = 1. - n*(2./nTheta_);
              scalar thetaBeg = Foam::acos(muBeg);
              scalar thetaEnd = Foam::acos(muEnd);
              deltaTheta = thetaEnd - thetaBeg;
              thetai = thetaBeg + (deltaTheta/2.); 
            }
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        scatter_, // ankur
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        // Currently 2D solution is limited to the x-y plane
        if (mesh_.solutionD()[vector::Z] != -1)
        {
            FatalErrorInFunction
                << "Currently 2D solution is limited to the x-y plane"
                << exit(FatalError);
        }

        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        nRay_ = 4*nPhi_;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi/(2.0*nPhi_);
        label i = 0;
        for (label m = 1; m <= 4*nPhi_; m++)
        {
            scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    absorptionEmission_,
                    scatter_, // ankur
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }
    // 1D
    else
    {
        // Currently 1D solution is limited to the x-direction
        if (mesh_.solutionD()[vector::X] != 1)
        {
            FatalErrorInFunction
                << "Currently 1D solution is limited to the x-direction"
                << exit(FatalError);
        }

        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        nRay_ = 2;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi;
        label i = 0;
        for (label m = 1; m <= 2; m++)
        {
            scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    absorptionEmission_,
                    scatter_, // ankur
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
    }

    // ankur, construct incident radiation for each band
    forAll(GLambda_, lambdaI)
    {
        GLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "GLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                G_
            )
        );
    }

    // ankur, construct energy fraction for each band
    forAll(enFracLambda_, lambdaI)
    {
        enFracLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "enFracLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("one", dimless, 1.0)
            )
        );
    }


    Info<< "fvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        Info<< '\t' << IRay_[rayId].I().name() << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << nl;
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::fvDOM(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qr", dimMass/pow3(dimTime), 0.0)
    ),
    Qem_
    (
        IOobject
        (
            "Qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qem", dimMass/pow3(dimTime), 0.0)
    ),
    Qin_
    (
        IOobject
        (
            "Qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimMass/pow3(dimTime), 0.0)
    ),
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
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    GLambda_(nLambda_), // ankur
    enFracLambda_(nLambda_), // ankur
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    uniSolidAngles_(coeffs_.lookupOrDefault<Switch>("uniformSolidAngles",false)) // ankur
{
    initialise();
}


Foam::radiation::fvDOM::fvDOM
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qr", dimMass/pow3(dimTime), 0.0)
    ),
    Qem_
    (
        IOobject
        (
            "Qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qem", dimMass/pow3(dimTime), 0.0)
    ),
    Qin_
    (
        IOobject
        (
            "Qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    GLambda_(nLambda_), // ankur
    enFracLambda_(nLambda_), // ankur
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    uniSolidAngles_(coeffs_.lookupOrDefault<Switch>("uniformSolidAngles",false)) // ankur
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::~fvDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOM::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::fvDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    absorptionEmission_->correctEnFrac(enFracLambda_, blackBody_); // ankur, to update energy fraction 
    
    // updateBlackBodyEmission();  // ankur, this not needed now, since the above call updates the fraction.. Also, the data for storing spectral fraction has now been defined in fvDOM class itself now.. 

    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0.0;
    label radIter = 0;
    do
    {
        Info<< "Radiation solver iter: " << radIter << endl;

        radIter++;
        maxResidual = 0.0;
        forAll(IRay_, rayI)
        {
            if (!rayIdConv[rayI])
            {
                scalar maxBandResidual = IRay_[rayI].correct();
                maxResidual = max(maxBandResidual, maxResidual);

                if (maxBandResidual < convergence_)
                {
                    rayIdConv[rayI] = true;
                }
            }
        }

    } while (maxResidual > convergence_ && radIter < maxIter_);

    updateG();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            // Only include continuous phase emission
            //4*absorptionEmission_->aCont()*physicoChemical::sigma  // ankur, commenting this line.. not correct for wide-band/WSGG models
            4*aDispCumm()*physicoChemical::sigma     //  ankur, aDispCumm function computes aDisp appropriately for wide-band/WSGG models
        )
    );
}


Foam::tmp<Foam::volScalarField::Internal> // ankur
Foam::radiation::fvDOM::Ru() const
{

//    const volScalarField::Internal& G =
//        G_();
//
//    const volScalarField::Internal E =
//        absorptionEmission_->ECont()()();
//
//    // Only include continuous phase absorption
//    const volScalarField::Internal a =
//        absorptionEmission_->aCont()()();
//
//    return a*G - E;

    // ankur
    tmp<volScalarField::Internal> tRu
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tRu",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    // ankur
    for (label j=0; j < nLambda_; j++)
    {
        tRu.ref() += absorptionEmission_->aCont(j)()()*GLambda_[j]() - absorptionEmission_->ECont(j)()();
    }

    return tRu; // ankur

}


void Foam::radiation::fvDOM::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}


// ankur
Foam::tmp<Foam::volScalarField>
Foam::radiation::fvDOM::aDispCumm() const
{
    tmp<volScalarField> aDispC
    (
        new volScalarField
        (
            IOobject
            (
                "aDispC",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    for (label j=0; j < nLambda_; j++)
    {
        aDispC.ref() += absorptionEmission_->aCont(j)*enFracLambda_[j]; 
    }   

    return aDispC;
}


void Foam::radiation::fvDOM::updateG()
{
    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    Qem_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
    Qin_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
    // ankur
    forAll(GLambda_,iLambda)
    {
    	GLambda_[iLambda] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    } 

    forAll(IRay_, rayI)
    {
        IRay_[rayI].addIntensity();
        G_ += IRay_[rayI].I()*IRay_[rayI].omega();
        qr_.boundaryFieldRef() += IRay_[rayI].qr().boundaryField();
        Qem_.boundaryFieldRef() += IRay_[rayI].Qem().boundaryField();
        Qin_.boundaryFieldRef() += IRay_[rayI].Qin().boundaryField();
        // ankur
        forAll(GLambda_,iLambda)
        {
          GLambda_[iLambda] += IRay_[rayI].ILambda(iLambda)*IRay_[rayI].omega();
        } 
    }
}


void Foam::radiation::fvDOM::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // assuming name is in the form: CHARS_rayId_lambdaId
    size_type i1 = name.find_first_of("_");
    size_type i2 = name.find_last_of("_");

    rayId = readLabel(IStringStream(name.substr(i1+1, i2-1))());
    lambdaId = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());
}

// ankur
Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::inScatEnergy(const label iLambda, const label sourDir) const
{
    tmp<volScalarField> tInScatEn
    (
        new volScalarField
        (
            IOobject
            (
                "inScatEn",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    forAll(IRay_, rayI)
    {
        tInScatEn.ref() += scatter_->pFunc(iLambda,rayI,sourDir)*IRay_[rayI].ILambda(iLambda)*IRay_[rayI].omega();

        //if (rayI != sourDir)
        //{
           //dimensionedScalar temp1("zero", dimLength, 1.0);
           //tInScatEn.ref() += scatter_->sigmaEff(iLambda)*IRay_[rayI].ILambda(iLambda)*IRay_[rayI].omega()/(4.*pi - IRay_[rayI].omega());
        //}
//Info << "scatIn: " << scatter_->pFunc(iLambda,rayI,sourDir).ref().primitiveFieldRef() << endl;
    }

    tInScatEn.ref() *= IRay_[sourDir].omega();

    return tInScatEn;
}

// ************************************************************************* //

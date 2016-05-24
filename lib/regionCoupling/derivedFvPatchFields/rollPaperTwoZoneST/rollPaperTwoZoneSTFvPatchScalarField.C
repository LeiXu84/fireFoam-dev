/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenCFD Ltd.
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

#include "rollPaperTwoZoneSTFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "turbulenceModel.H"
#include "thermoSingleLayerRollPaper.H"
#include "pyrolysisModel.H"

#include "constants.H"

#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const rollPaperTwoZoneSTFvPatchScalarField::filmModelType&
rollPaperTwoZoneSTFvPatchScalarField::
filmModel() const
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return *iter();
        }
    }


    FatalErrorIn
    (
        "const rollPaperTwoZoneSTFvPatchScalarField::"
        "filmModelType& "
        "rollPaperTwoZoneSTFvPatchScalarField::"
        "filmModel() const"
    )
        << "Unable to locate film region " << filmRegionName_
        << abort(FatalError);

    return **models.begin();
}


const rollPaperTwoZoneSTFvPatchScalarField::
pyrolysisModelType&
rollPaperTwoZoneSTFvPatchScalarField::
pyrModel() const
{
    HashTable<const pyrolysisModelType*> models =
        db().time().lookupClass<pyrolysisModelType>();

    forAllConstIter(HashTable<const pyrolysisModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == pyrolysisRegionName_)
        {
            return *iter();
        }
    }

    FatalErrorIn
    (
        "const rollPaperTwoZoneSTFvPatchScalarField::"
        "pyrolysisModelType& "
        "rollPaperTwoZoneSTFvPatchScalarField::"
        "pyrModel() const"
    )
        << "Unable to locate pyrolysis region " << pyrolysisRegionName_
        << abort(FatalError);

    return **models.begin();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rollPaperTwoZoneSTFvPatchScalarField::
rollPaperTwoZoneSTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    TnbrName_("undefined-Tnbr"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    neighbourFieldConvectiveName_("undefined-neigbourFieldConvectiveName"),
    fieldConvectiveName_("undefined-fieldConvectiveName"),
    KName_("undefined-K"),
    convectiveScaling_(1.0),
    convectiveCoefficient_(1.0),
    filmDeltaDry_(0.0),
    filmDeltaWet_(0.0),
    emissivity_(p.size(), 0.0),
    oldMode_(unknown)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


rollPaperTwoZoneSTFvPatchScalarField::
rollPaperTwoZoneSTFvPatchScalarField
(
    const rollPaperTwoZoneSTFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    radiationCoupledBase
    (
        p,
        psf.emissivityMethod(),
        psf.emissivity_
    ),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    fieldConvectiveName_(psf.fieldConvectiveName_),
    KName_(psf.KName_),
    convectiveScaling_(psf.convectiveScaling_),
    convectiveCoefficient_(psf.convectiveCoefficient_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    emissivity_(psf.emissivity_),
    oldMode_(psf.oldMode_)
{}


rollPaperTwoZoneSTFvPatchScalarField::
rollPaperTwoZoneSTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    radiationCoupledBase(p, dict),
    filmRegionName_
    (
        dict.lookupOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    pyrolysisRegionName_
    (
        dict.lookupOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    TnbrName_(dict.lookup("Tnbr")),
    neighbourFieldRadiativeName_(dict.lookup("neighbourFieldRadiativeName")),
    fieldRadiativeName_(dict.lookup("fieldRadiativeName")),
    neighbourFieldConvectiveName_(dict.lookup("neighbourFieldConvectiveName")),
    fieldConvectiveName_(dict.lookup("fieldConvectiveName")),
    KName_(dict.lookup("K")),
    convectiveScaling_(dict.lookupOrDefault<scalar>("convectiveScaling", 1.0)),
    convectiveCoefficient_(dict.lookupOrDefault<scalar>("convectiveCoefficient", 1.0)),
    filmDeltaDry_
    (
            dict.lookupOrDefault<scalar>("filmDeltaDry", 0.0000)
    ),
    filmDeltaWet_
    (
            dict.lookupOrDefault<scalar>("filmDeltaWet", 0.0002)
    ),
//    filmDeltaDry_(readScalar(dict.lookup("filmDeltaDry"))),
//    filmDeltaWet_(readScalar(dict.lookup("filmDeltaWet"))),
    emissivity_(p.size(), 0.0),
    oldMode_(unknown)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "rollPaperTwoZoneSTFvPatchScalarField::"
            "rollPaperTwoZoneSTFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


rollPaperTwoZoneSTFvPatchScalarField::
rollPaperTwoZoneSTFvPatchScalarField
(
    const rollPaperTwoZoneSTFvPatchScalarField&
        psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    radiationCoupledBase
    (
        psf.patch(),
        psf.emissivityMethod(),
        psf.emissivity_
    ),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    fieldConvectiveName_(psf.fieldConvectiveName_),
    KName_(psf.KName_),
    convectiveScaling_(psf.convectiveScaling_),
    convectiveCoefficient_(psf.convectiveCoefficient_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    emissivity_(psf.emissivity_),
    oldMode_(psf.oldMode_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rollPaperTwoZoneSTFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchI = patch().index();
    const label nbrPatchI = mpp.samplePolyPatch().index();
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchI];

    scalarField intFld(patchInternalField());

    const rollPaperTwoZoneSTFvPatchScalarField&
        nbrField =
        refCast
        <
            const rollPaperTwoZoneSTFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarField& Tp = *this;

    const scalarField K(this->kappa(*this));
    const scalarField nbrK(nbrField.kappa(*this));

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr(nbrK*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);

    scalarField nbrConvFlux = KDeltaNbr*(intFld - nbrIntFld);

    scalarField nbrTotalFlux = nbrConvFlux;
    // scalarList nbrRadField(nbrPatch.size(), 0.0);
    scalarList QrCoupled(nbrPatch.size(), 0.0);//kvm
    scalarList Tfilm(nbrPatch.size(), 0.0);//kvm
    scalarList alpha(nbrPatch.size(), 0.0);//kvm
    scalarList filmConv(nbrPatch.size(), 0.0);//kvm
    scalarList filmDelta(nbrPatch.size(), 0.0);//kvm
    scalarList myRadField(patch().size(), 0.0);
    
    scalarList convField(nbrPatch.size(), 0.0);//Ning
    scalarList TsurField(nbrPatch.size(), 0.0);//Ning

    const pyrolysisModelType& pyrolysis = pyrModel();
    const filmModelType& film = filmModel();

    // In solid
    if(neighbourFieldRadiativeName_ != "none") //nbr Radiation Qr
    {

        scalarField nbrConvFlux = convectiveCoefficient_*(intFld - nbrIntFld);
        
        const label filmPatchI =
            pyrolysis.nbrCoupledPatchID(film, patchI);

        const scalarField Qconvw(film.Qconvw(filmPatchI));

        // kvm, Qconvw is not right
        filmConv =
            // pyrRegion.mapRegionPatchField
            // (
            //     filmModel,
            //     patchI,
            //     filmPatchI,
            //     Qconvw,
            //     true
            // );
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "qFilmToWall",
                patchI,
                true
            );

        QrCoupled =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Qin",
                patchI,
                true
            );

        // Info << "QrCoupled " << QrCoupled << endl;

        Tfilm =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Tf",
                patchI,
                true
            );

        filmDelta =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "deltaf",
                patchI,
                true
            );
        alpha =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "alpha",
                patchI,
                true
            );

            const fvMesh& mesh = patch().boundaryMesh().mesh();
            if (! (mesh.foundObject<radiation::radiationModel>("radiationProperties")))
            {
                FatalErrorIn
                (
                    "rollPaperTwoZoneSTFvPatchScalarField::"
                    "rollPaperTwoZoneSTFvPatchScalarField\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const DimensionedField<scalar, volMesh>& iF,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "\n    radiationProperties file not found in pyrolysis region\n" 
                    << exit(FatalError);
            }
            const radiation::radiationModel& radiation =
                mesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );

            scalarField temissivity
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    //nbrPatch.index()
                    patch().index()
                ]
            );
            
            scalarField tabsorptivity
            (
                radiation.absorptionEmission().a()().boundaryField()
                [
                    //nbrPatch.index()
                    patch().index()
                ]
            );
            

        // nbrRadField =
        //     nbrPatch.lookupPatchField<volScalarField, scalar>
        //     (
        //         neighbourFieldRadiativeName_
        //     );

        convField =	//Ning
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldConvectiveName_
            );

        mpp.distribute(convField);	//Ning

        // Note: the Qr radiative flux is positive outgoing.
        // For a hot solid radiating into a cold fluid Qr will be negative.


        // Swap to obtain full local values of neighbour radiative heat flux
        // field
        // mpp.distribute(nbrRadField);

    //Qrad is negative going out of the solid
    //Qcon is positive going out of the solid

        // emissivity_ =
        //     patch().lookupPatchField<volScalarField, scalar>("emissivity");


        if (debug)
        {
            // scalar Qr = gSum(nbrRadField*patch().magSf());

            Info<< mesh.name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " :" << nl
                // << "     radiative heat  [W] : " << Qr << nl
//kvm                << "     predicted wallT [K] : " << gAverage(Twall) << nl
                << endl;
        }

        label nFixed = 0;

         // Estimate wetness of the film (1: wet , 0: dry)
         scalarField ratio
         (
            min
            (
                max
                (
                    (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
                    scalar(0.0)
                ),
                scalar(1.0)
            )
         );

	    //tmp<scalarField> kappaTmp = (kappa(*this));
	    //fvPatchField kappaTmp(kappa(*this));
	    
	    //-Ning: Add blockage effect
	    const volScalarField& blocking = db().lookupObject<volScalarField>("blockFactor");
	    const fvPatchScalarField blockingBoundary(blocking.boundaryField()[patch().index()]);
	    const scalarField blockingInternal = blockingBoundary.patchInternalField();
	    Info<<"Write net heat flux to the solid."<<endl;

	    const volScalarField& surfaceT = db().lookupObject<volScalarField>("Tsurface");
	    const fvPatchScalarField surfaceTBoundary(surfaceT.boundaryField()[patch().index()]);
	    const scalarField surfaceTInternal = surfaceTBoundary.patchInternalField();

	    const volScalarField& qheat = db().lookupObject<volScalarField>("Qnet");
	    fvPatchScalarField& qheatBoundary =
	    	const_cast<fvPatchScalarField&>(qheat.boundaryField()[patch().index()]);

        scalar qRadDryMax=-VGREAT;
        scalar qConvDryMax=-VGREAT;
        scalar qRadDryMin=VGREAT;
        scalar qConvDryMin=VGREAT;
        scalar qRadWetMax=-VGREAT;
        scalar qConvWetMax=-VGREAT;
        scalar qRadWetMin=VGREAT;
        scalar qConvWetMin=VGREAT;
        scalar qHeatBoundaryMax=-VGREAT;
        scalar qHeatBoundaryMin=VGREAT;

        const scalar sigma = constant::physicoChemical::sigma.value();

        //kvm-debug Info << "rollPaterTwoZoneST::qheatBoundary " << qheatBoundary << endl;

        forAll(*this, i)
        {
            scalar qConvWet = -filmConv[i];
            scalar qConvDry = -convField[i]; //Ning

            scalar qRadWet = 0.0; // all film absorption takes place in film model
            scalar qRadDry =
               -tabsorptivity[i]*QrCoupled[i]
               +temissivity[i]*sigma*pow4(surfaceTInternal[i]);	// count for emission from one side, if thin paper, need to add another side

            if(qRadDry > qRadDryMax)
            {
                qRadDryMax = qRadDry;
            }
            if(qConvDry > qConvDryMax)
            {
                qConvDryMax = qConvDry;
            }

            if(qRadDry < qRadDryMin)
            {
                qRadDryMin = qRadDry;
            }
            if(qConvDry < qConvDryMin)
            {
                qConvDryMin = qConvDry;
            }

            if(qRadWet > qRadWetMax)
            {
                qRadWetMax = qRadWet;
            }
            if(qConvWet > qConvWetMax)
            {
                qConvWetMax = qConvWet;
            }

            if(qRadWet < qRadWetMin)
            {
                qRadWetMin = qRadWet;
            }
            if(qConvWet < qConvWetMin)
            {
                qConvWetMin = qConvWet;
            }

            scalar qWet = qRadWet + qConvWet;
            scalar qDry = qRadDry + qConvDry;

            // kvm, this needs to account for wet surfaces, right?
            // alpha[i] = 0;
            qheatBoundary[i] = -((1.0-alpha[i])*(qDry) + alpha[i]*qWet);	    // Ning Two region Model
            //kvm-debug Info << "rollPaperTwoZoneST::qWet[" << i << "] " << qWet << endl;
            //kvm-debug Info << "rollPaperTwoZoneST::qDry[" << i << "] " << qDry << endl;
            //kvm-debug Info << "rollPaperTwoZoneST::alpha[" << i << "] " << alpha[i] << endl;

            if(qheatBoundary[i] < qHeatBoundaryMin)
            {
                qHeatBoundaryMin = qheatBoundary[i];
            }
            if(qheatBoundary[i] > qHeatBoundaryMax)
            {
                qHeatBoundaryMax = qheatBoundary[i];
            }


            scalar qDryBlockage = blockingInternal[i]*temissivity[i]*sigma
                *(pow4(operator[](i))-pow4(surfaceTInternal[i]))
                +(1.0-blockingInternal[i])*(qConvDry+qRadDry);

            nbrTotalFlux[i] = (1.0-alpha[i])*qDryBlockage + alpha[i]*qWet;
            //nbrTotalFlux[i] = 0.0;

            this->refValue()[i] = operator[](i);  // not used
            this->refGrad()[i] = -nbrTotalFlux[i]/K[i];                
            this->valueFraction()[i] = 0.0;
            nFixed++;
        }

        
        reduce(qRadDryMax,maxOp<scalar>());
        reduce(qConvDryMax,maxOp<scalar>());
        reduce(qRadDryMin,minOp<scalar>());
        reduce(qConvDryMin,minOp<scalar>());
        Info << "rollPaperTwoZoneST::qRadDryMin " << tab << db().time().timeName() << tab << qRadDryMin << endl;
        Info << "rollPaperTwoZoneST::qRadDryMax " << tab << db().time().timeName() << tab << qRadDryMax << endl;
        Info << "rollPaperTwoZoneST::qConvDryMin " << tab << db().time().timeName() << tab << qConvDryMin << endl;
        Info << "rollPaperTwoZoneST::qConvDryMax " << tab << db().time().timeName() << tab << qConvDryMax << endl;
        reduce(qRadWetMax,maxOp<scalar>());
        reduce(qConvWetMax,maxOp<scalar>());
        reduce(qRadWetMin,minOp<scalar>());
        reduce(qConvWetMin,minOp<scalar>());
        Info << "rollPaperTwoZoneST::qRadWetMin " << tab << db().time().timeName() << tab << qRadWetMin << endl;
        Info << "rollPaperTwoZoneST::qRadWetMax " << tab << db().time().timeName() << tab << qRadWetMax << endl;
        Info << "rollPaperTwoZoneST::qConvWetMin " << tab << db().time().timeName() << tab << qConvWetMin << endl;
        Info << "rollPaperTwoZoneST::qConvWetMax " << tab << db().time().timeName() << tab << qConvWetMax << endl;
        //debug Info << "QrCoupled " << max(QrCoupled) << endl;

        reduce(qHeatBoundaryMax,maxOp<scalar>());
        reduce(qHeatBoundaryMin,maxOp<scalar>());
        Info << "rollPaperTwoZoneST::qHeatBoundaryMin " << tab << db().time().timeName() << tab << qHeatBoundaryMin << endl;
        Info << "rollPaperTwoZoneST::qHeatBoundaryMax " << tab << db().time().timeName() << tab << qHeatBoundaryMax << endl;
        if (debug)
        {
            Pout<< "Using " << nFixed << " fixedValue out of " << this->size()
                << endl;
        }
    }
    else // In fluid
    {
//        typedef regionModels::pyrolysisModels::pyrolysisModel
//            pyrolysisModelType;
//
//        const regionModels::regionModel& pyrolysisRegion =
//            db().time().lookupObject<regionModels::regionModel>
//            (
//                "pyrolysisProperties"
//            );
//        const pyrolysisModelType& pyrolysisModel =
//            dynamic_cast<const pyrolysisModelType&>(pyrolysisRegion);
//
//        pyrolysisModelType& pyrolysis =
//            const_cast<pyrolysisModelType&>(pyrolysisModel);
//
//        const regionModels::regionModel& filmModel =
//            pyrolysisModel.db().time().lookupObject<regionModels::regionModel>
//            (
//                "surfaceFilmProperties"
//            );

        Tfilm =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Tf",
                nbrPatchI,
                true
            );

        filmDelta =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "deltaf",
                nbrPatchI,
                true
            );

        const scalarField& myRadField =
            patch().lookupPatchField<volScalarField, scalar>
            (
                fieldRadiativeName_
            );

        // do we still need to do mpp.distribute Tfilm
        mpp.distribute(Tfilm);

        mpp.distribute(filmDelta);

        TsurField =	//Ning
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                "Tsurface"
            );

        mpp.distribute(TsurField);	//Ning
        //kvm-debug Info << "rollPaperTwoZonesST::TsurField " << TsurField << endl;

        // use solid internal cell as Twall of gas. See if we can make re-radiation stable.
        // kvm, superseded by below  scalarField Twall = nbrIntFld;
        scalarList Twall(patch().size(), 0.0);//kvm
        
        // Estimate wetness of the film (1: wet , 0: dry)
        scalarField ratio
        (
           min
           (
               max
               (
                   (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
                   scalar(0.0)
               ),
               scalar(1.0)
           )
        );

        scalar TwallMax = -VGREAT;
        scalar TwallMin =  VGREAT;
        scalar TwallDryMin =  VGREAT;
        scalar TwallDryMax = -VGREAT;
        scalar TwallWetMin =  VGREAT;
        scalar TwallWetMax = -VGREAT;
        //kvm-debug Info << "rollPaperTwoZonesSt::Tfilm " << Tfilm << endl;
        forAll(*this, i)
        {
            scalar filmDeltaWet=0.0002; //m
            scalar filmDeltaDry=0.0000; //m

            scalar Twet = min(max(Tfilm[i], 298.15), 378.4);
            scalar Tdry = TsurField[i];
            // scalar Tdry = nbrIntFld[i];

            Twall[i] = ratio[i]*(Twet - Tdry) + Tdry;
            TwallWetMax = max(TwallWetMax,Twet);
            TwallWetMin = min(TwallWetMin,Twet);
            TwallDryMax = max(TwallDryMax,Tdry);
            TwallDryMin = min(TwallDryMin,Tdry);
            TwallMax = max(TwallMax,Twall[i]);
            TwallMin = min(TwallMin,Twall[i]);
        }
        reduce(TwallWetMax,maxOp<scalar>());
        reduce(TwallWetMin,minOp<scalar>());
        reduce(TwallDryMax,maxOp<scalar>());
        reduce(TwallDryMin,minOp<scalar>());
        reduce(TwallMax,maxOp<scalar>());
        reduce(TwallMin,minOp<scalar>());
        Info << "rollPaperTwoZonesSt::TwallWetMax " << tab << db().time().timeName() << tab << TwallWetMax << endl;
        Info << "rollPaperTwoZonesSt::TwallWetMin " << tab << db().time().timeName() << tab << TwallWetMin << endl;
        Info << "rollPaperTwoZonesSt::TwallDryMax " << tab << db().time().timeName() << tab << TwallDryMax << endl;
        Info << "rollPaperTwoZonesSt::TwallDryMin " << tab << db().time().timeName() << tab << TwallDryMin << endl;
        Info << "rollPaperTwoZonesSt::TwallMax " << tab << db().time().timeName() << tab << TwallMax << endl;
        Info << "rollPaperTwoZonesSt::TwallMin " << tab << db().time().timeName() << tab << TwallMin << endl;

        if (debug)
        {
            scalar Qr = gSum(myRadField*patch().magSf());

            // Info<< mesh.name() << ':'
            //     << patch().name() << ':'
            //     << this->internalField().name() << " :" << nl
            //     << "     radiative heat  [W] : " << Qr << nl
            //     << "     predicted wallT [K] : " << gAverage(Twall) << nl
            //     << endl;
        }

        this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0;
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        //scalar Qc = gSum(nbrConvFlux*patch().magSf());
        // kvm, this doesn't account for film convective...
        scalar Qc = gSum(convField*patch().magSf()); //Ning
        // scalar Qr = gSum(radField*patch().magSf());
        scalar Qt = gSum(nbrTotalFlux*patch().magSf());

        // Info<< mesh.name() << ':'
        //     << patch().name() << ':'
        //     << this->internalField().name() << " -> "
        //     << nbrMesh.name() << ':'
        //     << nbrPatch.name() << ':'
        //     << this->internalField().name() << " :"
        //     << " heatFlux:" << Qc
        //     // << " radiativeFlux:" << Qr
        //     << " totalFlux:" << Qt
        //     << " walltemperature "
        //     << " min:" << gMin(*this)
        //     << " max:" << gMax(*this)
        //     << " avg:" << gAverage(*this)
        //     << endl;
    }
}


void rollPaperTwoZoneSTFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>
    (
        os,
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    writeEntryIfDifferent<word>
    (
        os,
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldConvectiveName")<< //Ning
        neighbourFieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldConvectiveName")<< //Ning
        fieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldRadiativeName")<< //Ning
        fieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("K")<< //Ning
        KName_ << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
    radiationCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    rollPaperTwoZoneSTFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //

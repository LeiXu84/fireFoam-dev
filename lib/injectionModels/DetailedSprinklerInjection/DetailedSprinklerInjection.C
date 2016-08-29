/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "DetailedSprinklerInjection.H"
#include "TimeFunction1.H"
#include "mathematicalConstants.H"
#include "distributionModel.H"
#include "Pstream.H"
#include "OFstream.H"
#include "SortableList.H"

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

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DetailedSprinklerInjection<CloudType>::DetailedSprinklerInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
    )
    :
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    // operatingPressure_(readScalar(this->coeffDict().lookup("operatingPressure"))),
    radiusToSprinkler_(readScalar(this->coeffDict().lookup("radiusToSprinkler"))),
    positionList_(this->coeffDict().lookup("positionList")),
    nSprinklers_(positionList_.size()),
    nActivatedSprinklers_(0),
    injectorCellList_(nSprinklers_),
    tetFaceI_(-1),
    tetPtI_(-1),
    totalParcels_(0),
    direction_(this->coeffDict().lookup("direction")),
    centralCoreAngle_(this->coeffDict().template lookupOrDefault<scalar>("centralCoreAngle",90.0)),
    armDirection_(this->coeffDict().lookup("armDirection")),
    parcelDirVec_(direction_),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
        ),
    diameterCoefficient_(this->coeffDict().template lookupOrDefault<scalar>("diameterCoefficient",1.0)),
    momentumEfficiency_(this->coeffDict().template lookupOrDefault<scalar>("momentumEfficiency",0.8)),
    // flowRateProfile_
    // (
    //     TimeFunction1<scalar>
    //     (
    //         owner.db().time(),
    //         "flowRateProfile",
    //         this->coeffDict()

    //         )
    //     ),
    tanVec1_(vector::zero),
    tanVec2_(vector::zero),
    Tgas_
    (
        owner.db().objectRegistry::template lookupObject<volScalarField>
        (
            "T"
            )
        ),
    Ugas_
    (
        owner.db().objectRegistry::template lookupObject<volVectorField>
        (
            "U"
            )
        ),
    activeLinks_(this->coeffDict().subDict("rtiCoeffs"). template lookupOrDefault<Switch>("active",false)),
    RTI_(this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("RTI",200.0)),
    C_(this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("C",0.0)),
    initialTemperatureList_(nSprinklers_,298.15),
    activationTemperature_(this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("activationTemperature",432.0)),
    // have Andy implement this for Switch in SubModelBase.C
    // activated_(this->template getBaseProperty<Switch>("activated")),
    activatedList_(nSprinklers_,false),
    linkTemperatureList_(nSprinklers_,298.15),
    activationTimeList_(nSprinklers_,GREAT),
    rtiOutputInterval_(this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<label>("rtiOutputInterval",100)),
    filePtr_(),
    // lookup table related
    parcelIndex_(0),
    activeSprinklerIndex_(-1),
    injectionDeltaT_(0.),
    sampleSize_(0),
    tableName_("none"),
    sprinkler_("none"),
    orientation_("none"),
    nEle_(0),
    nAzi_(0),
    operatingPressure_(0.),
    tablePressures_(),
    kFactor_(0.),
    sampledRadius_(0.),
    idealFlowRate_(0.),
    rndGen_(0),
    cachedRndGen_(2,20000),
    phiTables_(),
    thetaTables_(),
    volFluxTables_(),
    dv50Tables_(),
    velMagTables_(),
    volFlux_(),
    volFlow_(),
    volFlowHist_(),
    diameterHist_(),
    particleHist_(),
    mappedIndices_(),
    dv50_(),
    velMag_(),
    area_(),
    avgVelMag_(),
    ele_(),
    azi_(),
    eleMin_(),
    aziMin_(),
    eleMax_(),
    aziMax_(),
    parcelParticles_(),
    sampleVolFlow_(),
    sampleD_(),
    sampleArea_(),
    sampleAvgVelMag_(),
    sampleEle_(),
    sampleAzi_(),
    sampleParcelParticles_()
{

    Info << "centralCoreAngle: " << centralCoreAngle_ << nl;

    readTableData();

    interpolatePressure();

    computeAreas();

    computeVolFlow();

    idealFlowRate_ = computeIdealFlowRate(); // L/s

    computeInjectionProperties();

    // scalar injectionDeltaT = 0.001;
    // sampleInjectionTable(injectionDeltaT);

    treatSprinklerActivation();

    // Normalise direction vector
    direction_ /= mag(direction_);

    // Dynamically set sampleSize_ based on time step
    // This will result in particles being injected at each timestep
    sampleSize_ = round(parcelsPerSecond_*injectionDeltaT_);
    /*DEBUG(injectionDeltaT_);*/
    /*DEBUG(parcelsPerSecond_);*/
    /*DEBUG(sampleSize_);*/
    setSampleSizes();

    totalParcels_ = sampleSize_;

    // writeVolumeFluxSprinklerInjectionProfile();

    tanVec1_ = armDirection_;
    tanVec2_ = direction_^tanVec1_;

    if(tanVec2_ == vector::zero)
    {
        FatalErrorIn
        (
            "DetailedSprinklerInjection()"
        )
            << "Sprinkler orientation " << direction_ << nl
            << "and arm orientation " << armDirection_ << nl
            << "are not consistent!" << nl
            << exit(FatalError);
    }

    // Set total volume to inject, gets overwritten later
    // this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_); // m3
    scalar conversionLiterPerM3 = 1000.0;
    this->volumeTotal_ = idealFlowRate_*conversionLiterPerM3*duration_; // m3
    /*DEBUG(this->volumeTotal_);*/

    cacheInjectorCells();

    if(activeLinks_){
        writeSprinklerHeader();
    }
}

template<class CloudType>
Foam::DetailedSprinklerInjection<CloudType>::DetailedSprinklerInjection
(
    const DetailedSprinklerInjection<CloudType>& im
    )
    :
    InjectionModel<CloudType>(im),
    duration_(im.duration_),
    // operatingPressure_(im.operatingPressure_),
    radiusToSprinkler_(im.radiusToSprinkler_),
    positionList_(im.positionList_),
    nSprinklers_(im.nSprinklers_),
    injectorCellList_(im.injectorCellList_),
    tetFaceI_(im.tetFaceI_),
    tetPtI_(im.tetPtI_),
    totalParcels_(im.totalParcels_),
    direction_(im.direction_),
    centralCoreAngle_(im.centralCoreAngle_),
    armDirection_(im.armDirection_),
    parcelDirVec_(im.parcelDirVec_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    diameterCoefficient_(im.diameterCoefficient_),
    momentumEfficiency_(im.momentumEfficiency_),
    // flowRateProfile_(im.flowRateProfile_),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec2_),
    Tgas_(im.Tgas_),
    Ugas_(im.Ugas_),
    activeLinks_(im.activeLinks_),
    RTI_(im.RTI_),
    C_(im.C_),
    initialTemperatureList_(im.initialTemperatureList_),
    activationTemperature_(im.activationTemperature_),
    activatedList_(im.activatedList_),
    linkTemperatureList_(im.linkTemperatureList_),
    activationTimeList_(im.activationTimeList_),
    filePtr_(im.filePtr_),
    // lookup table related
    parcelIndex_(im.parcelIndex_),
    activeSprinklerIndex_(im.activeSprinklerIndex_),
    injectionDeltaT_(im.injectionDeltaT_),
    sampleSize_(im.sampleSize_),
    tableName_(im.tableName_),
    sprinkler_(im.sprinkler_),
    orientation_(im.orientation_),
    nEle_(im.nEle_),
    nAzi_(im.nAzi_),
    operatingPressure_(im.operatingPressure_),
    tablePressures_(im.tablePressures_),
    kFactor_(im.kFactor_),
    sampledRadius_(im.sampledRadius_),
    idealFlowRate_(im.idealFlowRate_),
    rndGen_(im.rndGen_),
    cachedRndGen_(im.cachedRndGen_),
    phiTables_(im.phiTables_),
    thetaTables_(im.thetaTables_),
    volFluxTables_(im.volFluxTables_),
    dv50Tables_(im.dv50Tables_),
    velMagTables_(im.velMagTables_),
    volFlux_(im.volFlux_),
    volFlow_(im.volFlow_),
    volFlowHist_(im.volFlowHist_),
    diameterHist_(im.diameterHist_),
    particleHist_(im.particleHist_),
    mappedIndices_(im.mappedIndices_),
    dv50_(im.dv50_),
    velMag_(im.velMag_),
    area_(im.area_),
    avgVelMag_(im.avgVelMag_),
    ele_(im.ele_),
    azi_(im.azi_),
    eleMin_(im.eleMin_),
    aziMin_(im.aziMin_),
    eleMax_(im.eleMax_),
    aziMax_(im.aziMax_),
    parcelParticles_(im.parcelParticles_),
    sampleVolFlow_(im.sampleVolFlow_),
    sampleD_(im.sampleD_),
    sampleArea_(im.sampleArea_),
    sampleAvgVelMag_(im.sampleAvgVelMag_),
    sampleEle_(im.sampleEle_),
    sampleAzi_(im.sampleAzi_),
    sampleParcelParticles_(im.sampleParcelParticles_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DetailedSprinklerInjection<CloudType>::~DetailedSprinklerInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::DetailedSprinklerInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}

template<class CloudType>
Foam::scalar Foam::DetailedSprinklerInjection<CloudType>::timeStart()
{
    nActivatedSprinklers_=0;
    for(label i=0;i<nSprinklers_;i++)
    {
        if(activeLinks_)
        {
            if ( !activatedList_[i])
            {
                computeLinkTemperature(i);
            }
            else
            {
                nActivatedSprinklers_++;
            }
        }
        else
        {
            nActivatedSprinklers_++;
        }
    }

    for(label i=0;i<nSprinklers_;i++)
    {
        std::ostringstream buf;
        buf.precision(2);
        buf << i;

        if
            (
                this->owner().solution().transient()
                && this->owner().db().time().outputTime()
                )
        {
            this->template setBaseProperty<  bool  >
                (word("activated")+buf.str(),activatedList_[i]);
            this->template setBaseProperty< scalar >
                (word("activationTime")+buf.str(),activationTimeList_[i]);
            this->template setBaseProperty< scalar >
                (word("linkTemperature")+buf.str(),linkTemperatureList_[i]);
        }
    }

    if(activeLinks_){
        writeSprinklerData();
    }

    return this->SOI_;
}


template<class CloudType>
Foam::label Foam::DetailedSprinklerInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
    )
{

    // set injection time interval for later use
    // injectionDeltaT_ = this->injectionDeltaT_;
    injectionDeltaT_ = time1-time0;

    // Dynamically set sampleSize_ based on time step
    // This will result in particles being injected at each timestep
    sampleSize_ = round(parcelsPerSecond_*injectionDeltaT_);

    /*DEBUG(injectionDeltaT_);*/
    /*DEBUG(sampleSize_);*/
    setSampleSizes();
    totalParcels_ = sampleSize_;

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return totalParcels_*nActivatedSprinklers_ ;
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::DetailedSprinklerInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
    )
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        if(time1>time0){
            // return flowRateProfile_.integrate(time0, time1)*nActivatedSprinklers_;
            return this->volumeTotal_/duration_ * (time1 - time0) * nActivatedSprinklers_;
        }
    }
    else
    {
        return 0.0;
    }
    return 0.0;
}


template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::setParcelDirVec
(
    scalar elevationAngle,
    scalar azimuthalAngle
    )
{
    /*DEBUG("setParcelDirVec");*/
    const scalar deg2Rad = pi/180.0;
 
    // randomize quadrant
    label quadrant = -1;
    if(Pstream::master()){
        quadrant = rndGen_.integer(0,3);
    }
    reduce(quadrant,maxOp<scalar>());
    /*DEBUG(quadrant);*/
    switch(quadrant)
    {
        case 0:
            azimuthalAngle = azimuthalAngle;
            break;
        case 1:
            azimuthalAngle = 180.0-azimuthalAngle;
            break;
        case 2:
            azimuthalAngle = 180.0+azimuthalAngle;
            break;
        case 3:
            azimuthalAngle = 360.0-azimuthalAngle;
            break;
        default:
            DEBUG("something is wrong!");
    }
    /*DEBUG(azimuthalAngle);*/
    scalar alpha = cos(elevationAngle*deg2Rad);
    scalar dCorr = sin(elevationAngle*deg2Rad);
    vector normal = alpha*(tanVec1_*cos(azimuthalAngle*deg2Rad) + tanVec2_*sin(azimuthalAngle*deg2Rad));

    parcelDirVec_ = dCorr*direction_;
    parcelDirVec_ += normal;
    parcelDirVec_ /= mag(parcelDirVec_);

    return;
}


/* Called from InjectionModel::inject() in base class */
template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::setPositionAndCell
(
    const label parcelI_global,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
    )
{

    bool firstParcel=false;
    if(parcelI_global == 0){
        firstParcel = true;
    }

    if(firstParcel){
        // only sample table once per injection
        sampleInjectionTable(injectionDeltaT_);
    }

    // set local parcelI on a per sprinkler basis
    label whichActiveSprinkler=parcelI_global/totalParcels_;
    // If sprinklerIndex = 1, then the second active sprinkler should activate.
    // If sprinklerIndex = 2, then the third active sprinkler should activate.

    activeSprinklerIndex_=-1;
    label activeSprinklerCount=-1;
    for(label si=0;si<activatedList_.size();si++){
        if(activatedList_[si]){
            activeSprinklerCount++;
        }
        if(activeSprinklerCount==whichActiveSprinkler){
            activeSprinklerIndex_=si;
            break;
        }
    }

    label parcelI = parcelI_global%totalParcels_;
    parcelIndex_ = parcelI;

    // ideally, this should be a random value within the min/max of sample angle
    // scalar maxEleVar = nEle_/91.;
    // scalar maxAziVar = nAzi_/360.;
    scalar variationEle = -1.0;
    scalar variationAzi = -1.0;
    label index = mappedIndices_[parcelI];
    if(Pstream::master()){
        variationEle = (eleMax_[index]-eleMin_[index])*(rndGen_.scalar01()-0.0);
        variationAzi = (aziMax_[index]-aziMin_[index])*(rndGen_.scalar01()-0.0);
    }
    reduce(variationEle,maxOp<scalar>());
    reduce(variationAzi,maxOp<scalar>());
    // variationEle = 0.0;
    // variationAzi = 0.0;
    scalar elevationAngle = eleMin_[index]+variationEle; // sampleEle_[parcelI]+variationEle;
    scalar azimuthalAngle = aziMin_[index]+variationAzi; // sampleAzi_[parcelI]+variationAzi;
    
    /*Info << "DetailedInjection azi: " << aziMin_[index] << " " << azi_[index] << " " << aziMax_[index] << " " << variationAzi << " " << azimuthalAngle << endl;*/
    /*Allow central jet to cover a range of angles*/
    if (elevationAngle > centralCoreAngle_){
        elevationAngle = 90.0;
    }

    setParcelDirVec(elevationAngle,azimuthalAngle);

    position = positionList_[activeSprinklerIndex_]+radiusToSprinkler_*parcelDirVec_;

    this->findCellAtPosition
        (
         // the first four arguments are returned to InjectionModel::inject()
         cellOwner,
         tetFaceI,
         tetPtI,
         position,
         true
        );

    // DEBUG("exit setPositionAndCell");

    return;
}



template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::setParticleDiameter
(
    typename CloudType::parcelType& parcel
)
{

    // DEBUG("enter setParticleDiameter");

    parcel.d() = sampleD_[parcelIndex_];

    // DEBUG("exit setParticleDiameter");

    return;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::setParticleVelocity
(
    typename CloudType::parcelType& parcel
)
{
    // DEBUG("enter setParticleVelocity");

    parcel.U() = sampleAvgVelMag_[parcelIndex_] * parcelDirVec_;

    parcel.typeId() = activeSprinklerIndex_; // eventually gets overwritten

    // tmp<DimensionedField<scalar, volMesh> > sprinklerId (
    //     new DimensionedField<scalar, volMesh>
    //     (
    //         IOobject
    //         (
    //             this->name() + "sprinklerId" ,
    //             this->owner().time().timeName(),
    //             this->owner(),
    //             IOobject::READ_IF_PRESENT,
    //             IOobject::AUTO_WRITE
    //             ),
    //         this->owner().mesh(),
    //         dimensionedScalar("zero", dimMass, 0.0)
    //         )
    //     );
    // sprinklerId->write();


    // DEBUG("exit setParticleVelocity");

    return;
}

/* Called from InjectionModel::inject() in base class */
template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
    )
{
    // DEBUG("enter setProperties");

    // set particle diameter
    setParticleDiameter(parcel);

    setParticleVelocity(parcel);

    // DEBUG("exit setProperties");

    return;
}


template<class CloudType>
bool Foam::DetailedSprinklerInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::DetailedSprinklerInjection<CloudType>::validInjection(const label)
{
    return true;
}

template<class CloudType>
Foam::scalar Foam::DetailedSprinklerInjection<CloudType>::setNumberOfParticles
(
    const label parcels,
    const scalar volume,
    const scalar diameter,
    const scalar rho
    )
{

    scalar nP = sampleParcelParticles_[parcelIndex_];

    return nP;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::computeLinkTemperature
(
    const label sprinklerIndex
)
{

    scalar To = initialTemperatureList_[sprinklerIndex];
    // scalar Tg = 408.0;
    // scalar U = 1.0;
    // scalar dTg = 110.0;
    // scalar dt = 1.0;
    static List< scalar > dTeOld(nSprinklers_,0.0);
    static List< scalar > linkTemperature(nSprinklers_,To);


    // Initialize values on all procs to -GREAT;
    linkTemperature[sprinklerIndex]=-GREAT;
    scalar U=-GREAT;
    scalar Tg=-GREAT;

    const label ic = injectorCellList_[sprinklerIndex];
    if(ic > -1){
        U=mag(Ugas_[ic])+SMALL;
        Tg=Tgas_[ic];
        scalar dTg = Tg-To;
        const scalar deltaT = this->owner().time().deltaTValue();

        scalar dTe = sqrt(U)/RTI_*(dTg-(1+C_/(sqrt(U)+SMALL))*dTeOld[sprinklerIndex])*deltaT+dTeOld[sprinklerIndex];

        linkTemperature[sprinklerIndex]=initialTemperatureList_[sprinklerIndex]+dTe;


        dTeOld[sprinklerIndex]=dTe;
    }
    else{
        dTeOld[sprinklerIndex]=-GREAT;
    }
    // Reduce values to other nodes
    reduce(linkTemperature[sprinklerIndex], maxOp<scalar>());
    reduce(dTeOld[sprinklerIndex], maxOp<scalar>());
    reduce(U, maxOp<scalar>());
    reduce(Tg, maxOp<scalar>());

    linkTemperatureList_[sprinklerIndex] = linkTemperature[sprinklerIndex];
    if(linkTemperature[sprinklerIndex]>activationTemperature_){
        Info << "Sprinkler " << sprinklerIndex << " Activated!\n";
        activatedList_[sprinklerIndex] = true;
        activationTimeList_[sprinklerIndex] = this->owner().time().value();
        this->SOI_ = this->owner().time().value();
    }

    return;
}


template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::cacheInjectorCells()
{
    // Set/cache the injector cell
    // this is the location used for the rti calculation
    forAll(positionList_,index){

        this->findCellAtPosition
            (
                injectorCellList_[index],
                tetFaceI_,
                tetPtI_,
                positionList_[index],
                true
                );

    }

    return;
}


template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::writeVolumeFluxSprinklerInjectionProfile()
{
    // OFstream osP("sprinklerInjectionProfile");

    // osP << "#cell " << "\t";
    // osP << "e1 "    << "\t";
    // osP << "e2 "    << "\t";
    // osP << "a1 "    << "\t";
    // osP << "a2 "    << "\t";
    // osP << "area "  << "\t";
    // osP << "fr "    << "\t";
    // osP << "vfr "   << "\t";
    // osP << "npc "   << "\t";
    // osP << endl;

    // for(label iCell = 0; iCell < numberCells; iCell++){
    //     osP << iCell << "\t";
    //     osP << areaEachCell[iCell]  << "\t";
    //     osP << flowRateCell[iCell]*1000 << "\t"; // kg/s
    //     osP << volFCell_[iCell] << "\t";
    //     osP << endl;
    // }

    return;
}


template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::treatSprinklerActivation()
{
    if(activeLinks_){
        //sprinkler to be activated via RTI calculation
        this->SOI_=GREAT;
        scalar currentTime = this->owner().time().value();

        // read restart properties from <time>/uniform/lagrangian/reactingCloud1/reactingCloud1OutputProperties
        Switch atLeastOneActivated=false;
        for(label i=0;i<nSprinklers_;i++)
        {
            std::ostringstream buf;
            buf.precision(2);
            buf << i;

            activatedList_[i] = this->template getBaseProperty<  bool  >
                (word("activated")+buf.str());
            if(activatedList_[i] == true){
                atLeastOneActivated=true;
            }
            activationTimeList_[i] = this->template getBaseProperty< scalar >
                (word("activationTime")+buf.str());
            linkTemperatureList_[i] = this->template getBaseProperty< scalar >
                (word("linkTemperature")+buf.str());
            initialTemperatureList_[i]=linkTemperatureList_[i];
        }

        if(atLeastOneActivated && min(activationTimeList_) < currentTime)
        {
            this->SOI_ = min(activationTimeList_);
        }

        for(label i=0;i<nSprinklers_;i++){
            if(initialTemperatureList_[i]==0e0){
                // if linkTemperature not found in reactingCloud1OutputProperties then read from controlDict
                initialTemperatureList_[i] = this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("initialTemperature",298.0);
                activatedList_[i]=false;
                linkTemperatureList_[i] = initialTemperatureList_[i];
                activationTimeList_[i] = GREAT;
            }
        }
    }
    else{
        for(label i=0;i<nSprinklers_;i++){
            activatedList_[i] = true;
        }
    }

    return;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::writeSprinklerHeader()
{

    const fileName logDir = "sprinklerPostProcessing"/this->owner().time().timeName();
    mkDir(logDir);
    filePtr_ = new OFstream(logDir/"sprinklerRti.dat");
    OFstream& rtiData = *filePtr_;

    rtiData << "# number of sprinklers: " << nSprinklers_ << endl;
    rtiData << "# positions:" << endl;
    for(label i=0;i<nSprinklers_;i++){
        rtiData << "# sprinkler" << i << " ";
        rtiData << positionList_[i] << endl;
    }
    rtiData << "#t";
    for(label i=0;i<nSprinklers_;i++){
        rtiData << "     ";
        rtiData << "linkTemperature" << i;
    }
    rtiData << endl;

    return;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::writeSprinklerData()
{

    // write to file sprinklerRti.dat
    OFstream& rtiData = *filePtr_;

    static label counter = 0;
    counter++;
    if(counter%rtiOutputInterval_==0){
        rtiData.precision(5);
        rtiData.setf(std::ios::fixed);
        rtiData << this->owner().time().value();
        for(label i=0;i<nSprinklers_;i++){
            rtiData << "     ";
            rtiData << linkTemperatureList_[i];
        }
        rtiData << endl;
    }

    // write to log file
    Info << nl;
    Info << "Sprinkler injection:\n";
    Info << incrIndent;
    Info << indent << "Activated sprinklers: " << nActivatedSprinklers_ << " / " << nSprinklers_ << nl;
    Info << incrIndent;
    for(label i=0;i<nSprinklers_;i++){
        Info << indent << "sprinkler " << i << " (active, link temperature) : " << activatedList_[i] << " " << linkTemperatureList_[i] << nl;
    }
    Info << decrIndent << decrIndent;
    Info << nl;

    return;
}


//- Lookup table functions

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::readTableData()
{

    Info << "Reading sprinkler injection lookup table data\n";

    tableName_ = this->coeffDict().subDict("lookupTableCoeffs").template lookupOrDefault< word >("tableName","");
    fileName constant;
    if(Pstream::parRun()){
        constant =
            ".."/this->owner().db().time().constant();
    }
    else{
        Info << "tableName_: " << tableName_ << endl;
        constant =
            this->owner().db().time().constant();
    }

    operatingPressure_ = readScalar(this->coeffDict().subDict("lookupTableCoeffs").lookup("operatingPressure"));

    {
        IOdictionary dict
        (
            IOobject
            (
                tableName_,
                constant,
                this->owner().db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
        sprinkler_ = dict.template lookupOrDefault< word >("sprinkler","not defined");
        kFactor_ = readScalar(dict.lookup("kFactor")); // gpm/psi^0.5
        nEle_ = readLabel(dict.lookup("nEle"));
        nAzi_ = readLabel(dict.lookup("nAzi"));
        sampledRadius_ = readScalar(dict.lookup("radius"));
        Info << "sampledRadius_ " << sampledRadius_<< nl;
        orientation_ = dict.template lookupOrDefault< word >("orientation","not defined");

        /*Info << dict.toc();*/

        const dictionary& pressuresDict = dict.subDict("pressures");
        /*Info << pressuresDict.toc();*/
        const wordList pressureNames(pressuresDict.toc());
        /*Info << pressureNames << endl;*/
        forAll(pressureNames,i)
        {
            const dictionary& pressureSubDict = pressuresDict.subDict(pressureNames[i]);
            tablePressures_.append(readScalar(pressureSubDict.lookup("pressure")));
            /*Info << tablePressures_ << endl;*/

            phiTables_.append(pressureSubDict.lookup("phi"));
            thetaTables_.append(pressureSubDict.lookup("theta"));
            volFluxTables_.append(pressureSubDict.lookup("volFlux"));
            dv50Tables_.append(pressureSubDict.lookup("dv50"));
            velMagTables_.append(pressureSubDict.lookup("velMag"));

        }
        azi_.resize(phiTables_[0].size());
        ele_.resize(thetaTables_[0].size());
        azi_ = phiTables_[0]; // azi, phi
        ele_ = thetaTables_[0]; // ele, theta

    }

    return;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::interpolatePressure()
{

    Info << "Interpolate table data to operating pressure\n";

    // identify interpolation indices and weighting
    label begin = 0;
    label end = tablePressures_.size()-1;

    volFlux_.resize(volFluxTables_[begin].size());
    dv50_.resize(dv50Tables_[begin].size());
    velMag_.resize(velMagTables_[begin].size());

    if(tablePressures_.size()==1)
    {
        volFlux_ = volFluxTables_[begin];
        dv50_ = dv50Tables_[begin];
        velMag_ = velMagTables_[begin];
    }
    else if(operatingPressure_ <= tablePressures_[begin])
    {
        volFlux_ = (volFluxTables_[begin+1]-volFluxTables_[begin])/(tablePressures_[begin+1]-tablePressures_[begin])*(operatingPressure_ - tablePressures_[begin])+volFluxTables_[begin];
        dv50_ = (dv50Tables_[begin+1]-dv50Tables_[begin])/(tablePressures_[begin+1]-tablePressures_[begin])*(operatingPressure_ - tablePressures_[begin])+dv50Tables_[begin];
        velMag_ = (velMagTables_[begin+1]-velMagTables_[begin])/(tablePressures_[begin+1]-tablePressures_[begin])*(operatingPressure_ - tablePressures_[begin])+velMagTables_[begin];
    }
    else if(operatingPressure_ >= tablePressures_[end])
    {
        volFlux_ = (volFluxTables_[end]-volFluxTables_[end-1])/(tablePressures_[end]-tablePressures_[end-1])*(operatingPressure_ - tablePressures_[end])+volFluxTables_[end];
        dv50_ = (dv50Tables_[end]-dv50Tables_[end-1])/(tablePressures_[end]-tablePressures_[end-1])*(operatingPressure_ - tablePressures_[end])+dv50Tables_[end];
        velMag_ = (velMagTables_[end]-velMagTables_[end-1])/(tablePressures_[end]-tablePressures_[end-1])*(operatingPressure_ - tablePressures_[end])+velMagTables_[end];
    }
    else
    {
        for(label i=0;i<tablePressures_.size()-1;i++)
        {
            label i1 = i;
            label i2 = i+1;
            if(operatingPressure_ > tablePressures_[i1] && operatingPressure_ <= tablePressures_[i2])
            {
                scalar weight = (operatingPressure_-tablePressures_[i1])/(tablePressures_[i2]-tablePressures_[i1]);
                volFlux_ = (1.0-weight)*volFluxTables_[i1]+weight*volFluxTables_[i2];
                dv50_ = (1.0-weight)*dv50Tables_[i1]+weight*dv50Tables_[i2];
                velMag_ = (1.0-weight)*velMagTables_[i1]+weight*velMagTables_[i2];
            }
        }
    }

    /*Info << "volFlux: " << volFlux_ << endl;*/
    /*Info << "dv50: " << dv50_ << endl;*/
    /*Info << "velMag: " << velMag_ << endl;*/
    dv50_ = dv50_*1e-3; // convert dv50_ to m from mm

    /*DEBUG("leaving interpolatePressures()");*/
    return;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::computeAreas()
{
    // compute areas based on input radius
    
    /*DEBUG("computeAreas()");*/
    area_.resize(volFlux_.size());
    eleMin_.resize(ele_.size());
    eleMax_.resize(ele_.size());
    aziMin_.resize(azi_.size());
    aziMax_.resize(azi_.size());

    scalar totalArea = 0.0;
    
    scalar ele1;
    scalar ele2;
    scalar azi1;
    scalar azi2;
    for(label jazi=0; jazi < nAzi_; jazi++) // phi
    {
        for(label iele=0; iele < nEle_; iele++) // theta
        {
            /*DEBUG(jazi);*/
            /*DEBUG(iele);*/
            label iCell = jazi*(nEle_)+iele;
            /*DEBUG(iCell);*/
            if(iele==0)
            {
                // ele1 = ele_[iCell]; 
                ele1 = 0.0;
                ele2 = 0.5*(ele_[iCell]+ele_[iCell+1]); 
            }
            else if(iele==nEle_-1)
            {
                ele2 = ele_[iCell]; 
                ele1 = 0.5*(ele_[iCell]+ele_[iCell-1]); 
            }
            else
            {
                ele1 = 0.5*(ele_[iCell]+ele_[iCell-1]); 
                ele2 = 0.5*(ele_[iCell]+ele_[iCell+1]); 
            }
            if(jazi==0)
            {
                azi1 = azi_[iCell]; 
                azi2 = 0.5*(azi_[iCell]+azi_[iCell+1*nEle_]); 
            }
            else if(jazi==nAzi_-1)
            {
                azi2 = azi_[iCell]; 
                azi1 = 0.5*(azi_[iCell]+azi_[iCell-1*nEle_]); 
            }
            else
            {
                azi1 = 0.5*(azi_[iCell]+azi_[iCell-1*nEle_]); 
                azi2 = 0.5*(azi_[iCell]+azi_[iCell+1*nEle_]); 
            }
            /*DEBUG(ele_[iCell]);*/
            /*DEBUG(azi_[iCell]);*/
            /*DEBUG(ele1);*/
            /*DEBUG(ele2);*/
            /*DEBUG(azi1);*/
            /*DEBUG(azi2);*/
            eleMin_[iCell] = ele1;
            eleMax_[iCell] = ele2;
            aziMin_[iCell] = azi1;
            aziMax_[iCell] = azi2;

            const scalar deg2Rad = mathematical::pi/180.0;

            area_[iCell] = sampledRadius_*sampledRadius_*
                            (sin(deg2Rad*ele2)-sin(deg2Rad*ele1))*(deg2Rad*azi2-deg2Rad*azi1);

            totalArea = totalArea+area_[iCell];
        }
    }

    /*DEBUG(totalArea);*/
    scalar idealArea = 4.*mathematical::pi*pow(sampledRadius_,2)/2./4.;
    /*DEBUG(idealArea);*/

    /*DEBUG("leaving computeAreas()");*/
    return;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::computeVolFlow()
{

    /*DEBUG("computeVolFlow()");*/
    /*Info << volFlux_ << endl;*/
    // volFlow_.resize(volFlux_.size());
    volFlow_ = volFlux_*area_;
    volFlowHist_.resize(volFlow_.size(),0.0);
    diameterHist_.resize(volFlow_.size());
    particleHist_.resize(volFlow_.size());
    
    /*Info << "volFlow_: " << volFlow_ << endl;*/
    /*Info << "volFlow_.indices(): " << volFlow_.indices() << endl;*/
    // sort volFlow_
    volFlow_.sort();
    /*Info << "volFlow_: " << volFlow_ << endl;*/
    /*Info << "volFlow_.indices(): " << volFlow_.indices() << endl;*/
    /*DEBUG("leaving computeVolFlow()");*/

}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::computeInjectionProperties()
{

    // compute sprinker orifice diameter
    // K = 29.83 c D^2
    // where D is in inches, K in british units
    scalar Cd = 0.93; // tapered orifice, orifice coeficient
    scalar conversionKFactor = 14.464; //  gpm/psi^0.5 per lpm/bar^0.5
    scalar conversionMeterPerInch = 0.0254; //  meter per inch
    scalar conversionPsiPerBar = 14.503774; // psi per bar
    scalar factor = 29.83/pow(conversionMeterPerInch,2)*conversionKFactor; // L/min per m2
    Info << "factor: " << factor << endl;
    scalar orificeDiameter = sqrt(kFactor_*conversionKFactor/factor/Cd); // m
    Info << "orificeDiameter: " << orificeDiameter << endl;

    // compute velocity based on U = mdot/rho/A
    scalar rhow = 1000.0; // density of water, kg/m3
    //old scalar uJet = kFactor_*sqrt(operatingPressure_)/15.0/rhow/pi/pow(orificeDiameter,2.0); // m
    scalar Q = kFactor_*conversionKFactor*sqrt(operatingPressure_/conversionPsiPerBar); // lpm
    Info << "Q: " << Q << endl;
    Info << "operatingPressure_ (psig): " << operatingPressure_ << endl;
    Info << "kFactor_ (gpm/psi^0.5): " << kFactor_ << endl;

    scalar mdot = Q/60.;  // kg/s (assuming 1 l = 1 kg)
    scalar area = pi*pow(orificeDiameter,2.0)/4.0; // m2
    scalar uJet = mdot/(rhow*area); // m/s
    Info << "uJet: " << uJet << endl;
    Info << "diameterCoefficient: " << diameterCoefficient_ << endl;
    Info << "momentumEfficiency: " << momentumEfficiency_ << endl;
    // velMag_ = uJet*momentumEfficiency_; // account for momentum loss during atomization process

    // compute dv50
    scalar C = 1.9; // sprinkler specific empirical factor for determining dv50
    scalar sigmaw = 72.8e-3; // N/m @ 20 deg C
    scalar We = rhow*uJet*uJet*orificeDiameter/sigmaw;
    Info << "We: " << We << endl;

    // dv50_ = C*orificeDiameter/pow(We,0.33333); // m
    // dv50_ = diameterCoefficient_*dv50_;
    // Info << "dv50: " << dv50_ << endl;

    // compute particles per parcel
    // I don't think the following is needed
    parcelParticles_.resize(azi_.size());
    forAll(parcelParticles_,nc){
        scalar volumetricFlowRate = volFlux_[nc]*area_[nc]; // L/s
        scalar injectionDeltaT = 1.0; // s, representative only
        scalar conversionLiterPerM3 = 1000.0;
        scalar volumeToInject = volumetricFlowRate*injectionDeltaT/conversionLiterPerM3; // m3
        // scalar radius = 0.5*dv50_[nc]; // m
        scalar radius = 0.5*dv50_[0]; // m
        scalar singleVolume = 4./3.*pi*pow(radius,3); // m3
        parcelParticles_[nc] = volumeToInject/singleVolume; // None
    }

    return;
}

template<class CloudType>
Foam::label Foam::DetailedSprinklerInjection<CloudType>::weightedSampling
()
{
    // weighted sampling
    scalar flowSum = 0.0;
    flowSum = sum(volFlow_);

    scalar y(-1.);
    if(Pstream::master()){
        y = rndGen_.scalar01()*flowSum;
    }
    reduce(y,maxOp<scalar>());
    // DEBUG(y);

    forAll(volFlow_,i)
    {
        if(y < volFlow_[i])
            return i;
        y -= volFlow_[i];
    }
    // static_assert(false,"should never get here");

}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::sampleInjectionTable
(
    const scalar injectionDeltaT
)
{

    /*DEBUG("here");*/

    // sample (sampleSize_) number of indices
    mappedIndices_.resize(sampleSize_);
    labelList histogram(volFlow_.size(),0);
    forAll(mappedIndices_,j)
    {
        label i = weightedSampling();
        /*DEBUG(i);*/
        histogram[i]++;
        mappedIndices_[j] = volFlow_.indices()[i];
    }
    /*Info << "histogram: " << histogram << endl;*/
    /*Info << "mappedIndices: " << mappedIndices_ << endl;*/

   // knowing sampleSize_ and idealFlowRate_ then compute equal
   // flow for each package and assign to corresponding mappedIndices

    scalar parcelFlow = idealFlowRate_/sampleSize_; // L/s per parcel
     
    /*DEBUG(parcelFlow);*/
    /*DEBUG(idealFlowRate_);*/
    
    labelList sampleIndices(mappedIndices_);


    scalar sumVolumetricFlowRate(0.0);

    // needed if procs out of sync with random number generation
    // if(Pstream::master()){
    // }
    // else{
    //     forAll(sampleIndices,i){
    //         sampleIndices[i] = -1;
    //     }
    // }

    // forAll(sampleIndices,i){
    //     reduce(sampleIndices[i],maxOp<label>());
    // }

    for(int nc=0;nc<sampleSize_;nc++){
        label index = sampleIndices[nc];
        sampleAzi_[nc] = azi_[index];
        sampleEle_[nc] = ele_[index];
        sampleVolFlow_[nc] = parcelFlow; // volFlux_[index];
        scalar dt = this->owner().time().deltaT().value();
        volFlowHist_[index] += parcelFlow*dt;
        sampleArea_[nc] = area_[index];
        // sampleD_[nc] = dv50_[linearIndex]; // m
        // sampleD_[nc] = dv50_; // m
        sampleD_[nc] = sampleRosinRammler(dv50_[index]); // m
        diameterHist_[index].append(sampleD_[nc]);
        // sampleAvgVelMag_[nc] = avgVelMag_[linearIndex]; // m/s
        sampleAvgVelMag_[nc] = momentumEfficiency_*velMag_[index]; // m/s

        // scalar volumetricFlowRate = sampleVolFlow_[nc]*sampleArea_[nc];
        sumVolumetricFlowRate += sampleVolFlow_[nc]; // L/s
    }

    scalar ratio = sumVolumetricFlowRate/idealFlowRate_;

    // sumVolumetricFlowRate = 0.0;
    // scalar sumParcelParticles = 0.0;

    // forAll(volFlowHist_,i)
    // {
    //     Info << "azi: " << azi_[i] << " ele: " << ele_[i] << " volFlowHist: " << volFlowHist_[i] << " volFlow: " << volFlux_[i]*area_[i] << endl;
    // }

    scalar conversionLiterPerM3 = 1000.0;
    forAll(sampleIndices,nc){
        label index = sampleIndices[nc];
        // scalar volumetricFlowRate = sampleVolFlow_[nc]*sampleArea_[nc]/ratio; // L/s
        // new approach doesn't need to normalize
        // scalar volumetricFlowRate = sampleVolFlow_[nc]/ratio; // L/s
        // sumVolumetricFlowRate += volumetricFlowRate; // L/s

        scalar volumeToInject = sampleVolFlow_[nc]*injectionDeltaT/conversionLiterPerM3; // m3
        scalar radius = 0.5*sampleD_[nc]; // m
        scalar singleVolume = 4./3.*pi*pow(radius,3); // m3
        sampleParcelParticles_[nc] = volumeToInject/singleVolume; // None
        particleHist_[index].append(sampleParcelParticles_[nc]);
        // sumParcelParticles += sampleParcelParticles_[nc];
    }

    // forAll(diameterHist_,i)
    // {
    //     Info << "azi: " << azi_[i] << " ele: " << ele_[i] << " volFlowHist: " << volFlowHist_[i] << " volFlow: " << volFlux_[i]*area_[i] << endl;
    //     Info << "dv50: " << dv50_[i] << endl << "diameterHist: " << diameterHist_[i] << " particleHist: " << particleHist_[i] << endl;
    //     Info << "sampled-dv50: " << computeDv50(diameterHist_[i],particleHist_[i]) << " " << dv50_[i] << endl;
    // }

    return;
}


template<class CloudType>
Foam::scalar Foam::DetailedSprinklerInjection<CloudType>::computeIdealFlowRate()
{
    scalar secondsPerMinute = 60.0;

    scalar conversionKFactor = 14.464; //  gpm/psi^0.5 per lpm/bar^0.5
    scalar conversionPsiPerBar = 14.503774; // psi per bar
    scalar gpm_per_lps = 15.85; // gpm per l/s

    scalar idealFlowRate = kFactor_*conversionKFactor*sqrt(operatingPressure_/conversionPsiPerBar); // L/min
    idealFlowRate /= secondsPerMinute; // L/s
    Info << "idealFlowRate " << (idealFlowRate) << " L/s " << nl;
    Info << "idealFlowRate " << (idealFlowRate*gpm_per_lps) << " gpm " << nl;

    return idealFlowRate;
}

template<class CloudType>
void Foam::DetailedSprinklerInjection<CloudType>::setSampleSizes()
{

    sampleVolFlow_.resize(sampleSize_);
    sampleD_.resize(sampleSize_);
    sampleArea_.resize(sampleSize_);
    sampleAvgVelMag_.resize(sampleSize_);
    sampleEle_.resize(sampleSize_);
    sampleAzi_.resize(sampleSize_);
    sampleParcelParticles_.resize(sampleSize_);

    return;
}

template<class CloudType>
Foam::scalar Foam::DetailedSprinklerInjection<CloudType>::sampleRosinRammler(const scalar dv50)
{
    scalar n_ = 2.6; // as recommeded by fds :)
    scalar d_ = dv50/pow(0.693,1./n_); // d_ is D_{0.632}
    /*Info << "d_ " << d_ << endl;*/
    scalar maxValue_ = 1.0*d_*pow(6.9077,1./n_); // D_{0.999}
    /*Info << "maxValue_ " << maxValue_ << endl;*/
    // minValue_ must be small enough (0.001*D_{0.1}) in order for sampling to return dv50
    scalar minValue_ = 0.001*d_*pow(0.1054,1./n_); // 0.1*D_{0.1}
    /*Info << "minValue_ " << minValue_ << endl;*/
    scalar K = 1.0 - exp(-pow((maxValue_ - minValue_)/d_, n_));
    // why cached?
    scalar y = cachedRndGen_.sample01<scalar>();
    // scalar y = rndGen_.scalar01();
    scalar x = minValue_ + d_*::pow(-log(1.0 - y*K), 1.0/n_);
    return x;
}

template<class CloudType>
Foam::scalar Foam::DetailedSprinklerInjection<CloudType>::computeDv50(const scalarList& d,const scalarList& np)
{
    SortableList< scalar > dsort;
    dsort = 1.0*d;
    dsort.sort();
    scalarList r = dsort/2.;
    scalarList volume = 4./3.*pi*pow(r,3); // m3
    scalarList cumVol(volume.size(),0.0);
    labelList mapping = dsort.indices();
    forAll(volume,i)
    {
        volume[i] *= np[mapping[i]];
        Info << "volume[" << i << "]: " << volume[i] << endl;
    }
    forAll(cumVol,i)
    {
        if(i==0)
        {
            cumVol[i] = volume[i];
        }
        else
        {
            cumVol[i] += volume[i]+cumVol[i-1];
        }
    }
    Info << cumVol << endl;
    scalar v50 = 0.5*sum(volume);
    scalar x(0.0);    
    forAll(dsort,i)
    {
        if(cumVol[i]>v50)
        {
            x = dsort[i];
            break;
        }
    }
    return x;
}

// ************************************************************************* //

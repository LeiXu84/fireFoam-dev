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

  Class
  Foam::DetailedSprinklerInjection2

  Description
  Detailed-Sprinkler injection

  - User specifies
  - time of start of injection
  - injector positions
  - direction (along injection axis)
  - parcel flow rate
  - parcel velocity

  SourceFiles
  DetailedSprinklerInjection2.C

  \*---------------------------------------------------------------------------*/

#ifndef DetailedSprinklerInjection2_H
#define DetailedSprinklerInjection2_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "Random.H"
#include "TimeFunction1.H"
#include "OFstream.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

    template<class Type>
    class TimeFunction1;

    class distributionModel;

/*---------------------------------------------------------------------------*\
  Class DetailedSprinklerInjection2 Declaration
  \*---------------------------------------------------------------------------*/

    template<class CloudType>
    class DetailedSprinklerInjection2
        :
        public InjectionModel<CloudType>
    {
    private:

        // Private data

        //- Injection duration [s]
        const scalar duration_;
        //- Operating pressure [bar]
        // const scalar operatingPressure_;

        const scalar radiusToSprinkler_;

        //- List of injector positions [m]
        List< vector > positionList_;

        //- Total number of sprinklers
        label nSprinklers_;

        //- Total number of activated sprinklers
        label nActivatedSprinklers_;

        //- List of cells containing injector positions []
        List< label > injectorCellList_;

        //- Index of tet face for injector cell
        label tetFaceI_;

        //- Index of tet point for injector cell
        label tetPtI_;

        //- total number of injected parcels in each time step
        label totalParcels_;

        //- Sprinkler Spray Inject Center direction []
        vector direction_;
        
        //- Minimum angle of central core jet
        scalar centralCoreAngle_;

        //- One Sprinkler Frame Arm direction []
        vector armDirection_;
        vector parcelDirVec_;
        vector parcelDirVecVelocity_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Coefficient to scale the dv50 diameter []
        scalar diameterCoefficient_;

        //- Momentum efficiency []
        scalar momentumEfficiency_;

        //- Elevation angle deviation [deg]
        scalar elevationAngleDeviation_;

        //- Volume flow reduction factor []
        scalar reductionFactor_;

        //- Flow rate profile relative to SOI []
        // const TimeFunction1<scalar> flowRateProfile_;

        // Tangential vectors to the direction vector

        //- First tangential vector
        vector tanVec1_;

        //- Second tangential vector
        vector tanVec2_;

        const volScalarField& Tgas_;
        const volVectorField& Ugas_;
        Switch activeLinks_;
        scalar RTI_;
        scalar RTI_deflector_;
        scalar C_;
        List< scalar > initialTemperatureList_;
        scalar activationTemperature_;
        List< Switch > activatedList_;
        List< label >  activatedLabelList_;
        List< scalar > linkTemperatureList_;
        List< scalar > activationTimeList_;
        label rtiOutputInterval_;

        //- pointer to log file for sprinkler diagnostic output
        OFstream *filePtr_;


        //- lookup table related
        Switch verification_;
        Switch bernoulli_;
        label parcelIndex_;
        label activeSprinklerIndex_;
        scalar injectionDeltaT_;
        label sampleSize_;
        word tableName_;
        word sprinkler_;
        word orientation_;
        List< scalar > phiMinMax_;
        List< scalar > thetaMinMax_;
        label nEle_;
        label nAzi_;
        scalar operatingPressure_;
        scalar uJet_;
        scalarList tablePressures_;
        scalarList velocityCorrections_;
        scalar  velocityCorrection_;
        //K-factor of operating pressure lpm/bar^0.5
        scalar kFactor_;
        scalar ratio_;
        scalar totalTime_;
        scalar sampledRadius_;
        scalar idealFlowRate_;
        Random rndGen_;
        cachedRandom cachedRndGen_;

        // Table values
        List< List< scalar > > volFluxTables_;
        List< List< scalar > > dv50Tables_;
        List< List< scalar > > velMagTables_;

        // Interpolated values
        List< scalar > volFlux_;
        SortableList< scalar > volFlow_;
        List< scalar > volFlowHist_;
        List< List< scalar > > diameterHist_;
        List< List< scalar > > particleHist_;
        labelList mappedIndices_;
        List< scalar > dv50_;
        List< scalar > velMag_;
        // scalar dv50_;
        // scalar velMag_;
        List< scalar > area_;
        List< scalar > avgVelMag_;
        List< scalar > ele_;
        List< scalar > azi_;
        List< scalar > eleMin_;
        List< scalar > aziMin_;
        List< scalar > eleMax_;
        List< scalar > aziMax_;
        List< scalar > parcelParticles_;
        
        List< scalar > sampleVolFlow_;
        List< scalar > sampleD_;
        List< scalar > sampleArea_;
        List< scalar > sampleAvgVelMag_;
        List< scalar > sampleEle_;
        List< scalar > sampleAzi_;
        List< scalar > sampleParcelParticles_;


    protected:

        // sets the diameter of each injected parcel
        void setParticleDiameter(    
            typename CloudType::parcelType& parcel
            );

        // sets the velocity of each injected parcel
        void setParticleVelocity(    
            typename CloudType::parcelType& parcel
            );

        scalar getElevationAngle(Random& rndGen);
        scalar getAzimuthalAngle(Random& rndGen);
        void   setParcelDirVec(scalar elevationAngle,scalar azimuthalAngle);


    public:

        //- Runtime type information
        TypeName("detailedSprinklerInjection2");


        // Constructors

        //- Construct from dictionary
        DetailedSprinklerInjection2
        (
            const dictionary& dict, 
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        DetailedSprinklerInjection2(const DetailedSprinklerInjection2<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
                (
                    new DetailedSprinklerInjection2<CloudType>(*this)
                );
        }


        //- Destructor
        virtual ~DetailedSprinklerInjection2();


        // Member Functions

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Return the start-of-injection time
        virtual scalar timeStart();

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

        //- Set the injection position and owner cell
        virtual void setPositionAndCell
        (
            const label parcelI,
            const label nParcels,
            const scalar time,
            vector& position,
            label& cellOwner,
            label& tetFaceI,
            label& tetPtI
            );

            
        //- Set the parcel properties
        virtual void setProperties
        (
            const label parcelI,
            const label nParcels,
            const scalar time,
            typename CloudType::parcelType& parcel
            );

        //- Flag to identify whether model fully describes the parcel
        virtual bool fullyDescribed() const;

        //- Return flag to identify whether or not injection of parcelI is
        //  permitted
        virtual bool validInjection(const label parcelI);


        //- Set number of particles to inject given parcel properties
        virtual scalar setNumberOfParticles
        (
            const label parcels,
            const scalar volume,
            const scalar diameter,
            const scalar rho
            );

        void computeLinkTemperature(const label injectorCell);
        
        void cacheInjectorCells();

        void writeVolumeFluxSprinklerInjectionProfile();

        void treatSprinklerActivation();

        void writeSprinklerHeader();

        void writeSprinklerData();

        //- Table-lookup functions
        void readTableData();
        void interpolatePressure();
        void computeAreas();
        void computeVolFlow();
        void computeInjectionProperties();
        label weightedSampling();
        void sampleInjectionTable(const scalar injectionDeltaT);
        scalar computeIdealFlowRate();
        void setSampleSizes();
        scalar sampleRosinRammler(const scalar dv50);
        scalar computeDv50(const scalarList& d,const scalarList& np);


    // Public data

//        //- Phase type enumeration
//        enum lookupDataType
//        {
//            VELOCITY,
//            DIAMETER,
//            VOLUMETRIC_FLUX,
//            UNKNOWN
//        };
//
//        //- Corresponding word representations for lookup data type enumerations
//        static const NamedEnum<lookupDataType, 4> lookupDataTypeNames_;
//
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DetailedSprinklerInjection2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

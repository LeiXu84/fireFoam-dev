/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::ThermoCloud

Description
    Templated base class for thermodynamic cloud

    - Adds to kinematic cloud
      - Heat transfer

SourceFiles
    ThermoCloudI.H
    ThermoCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoCloud_H
#define ThermoCloud_H

#include "KinematicCloud.H"
#include "thermoCloud.H"
#include "SLGThermo.H"
#include "radiationModel.H" // ankur

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class HeatTransferModel;

/*---------------------------------------------------------------------------*\
                        Class ThermoCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermoCloud
:
    public CloudType,
    public thermoCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef ThermoCloud<CloudType> thermoCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<ThermoCloud<CloudType>> cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        ThermoCloud(const ThermoCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ThermoCloud&);


protected:

    // Protected data

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the carrier gas fields

            //- SLG thermodynamics package
            const SLGThermo& thermo_;

            //- Temperature [K]
            const volScalarField& T_;

            //- Pressure [Pa]
            const volScalarField& p_;


        // References to the cloud sub-models

            //- Heat transfer model
            autoPtr<HeatTransferModel<ThermoCloud<CloudType>>>
                heatTransferModel_;


        // Reference to the particle integration schemes

            //- Temperature integration
            autoPtr<scalarIntegrationScheme> TIntegrator_;


        // Modelling options

            //- Include radiation
            Switch radiation_;

            //- Include coupled radiation
            Switch coupledRadiation_; // ankur

            //- Reference to radiation model
            //const radiation::radiationModel& radModel_;

            //- Number of wavelengths/bands
            label nBands_; // ankur

            //- Radiation property type
            word radProp_; // ankur

            //- Constant radiation properties
            scalar constAbsEff_; // ankur
            scalar constSctEff_; // ankur

            //- Diameter dependent, banded radiation properties
            label numPropBands_; // ankur
            scalar beamLen_; // ankur
            List<scalar> energyFrac_; // ankur
            List<scalar> diaVal_; // ankur
            List<List<scalar>> absEff_; // ankur
            List<List<scalar>> sctEff_; // ankur
            List<List<scalar>> asyFac_; // ankur
 
            //- Radiation sum of parcel projected areas
            PtrList<volScalarField::Internal> radAreaP_; // ankur

            //- Radiation sum of parcel projected areas for scattering
            PtrList<volScalarField::Internal> radAreaPSc_; // ankur

            //- Radiation sum of parcel temperature^4
            PtrList<volScalarField::Internal> radT4_; // ankur

            //- Radiation sum of parcel projected areas * temperature^4
            PtrList<volScalarField::Internal> radAreaPT4_; // ankur

            //- Radiation sum of scattering coefficient * asymmetry factor
            PtrList<volScalarField::Internal> radAreaPScAsy_; // ankur

            //- Normalization factor for in-scattering terms  (as a function of bands and space)
            //- Assuming that this factor is same for different directions (different rays) for now..
            PtrList<volScalarField::Internal> CScat_; // ankur
   
            //- Angular-discretization dependent constants appearing in the expression for computation of CScat
            //- Computed only once at the start of the run
            List<scalar> CScatCoeffs_; // ankur

            //- Bool for setting CScat
            bool setCScat_;  // ankur

            //- Bool for setting CScatCoeffs
            bool setCScatCoeffs_;  // ankur

        // Sources

            //- Sensible enthalpy transfer [J/kg]
            autoPtr<volScalarField::Internal> hsTrans_;

            //- Coefficient for carrier phase hs equation [W/K]
            autoPtr<volScalarField::Internal> hsCoeff_;


    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();

            //-Set CScatCoeffs
            void setCScatCoeffs(); // ankur

            //-Set CScat
            void setCScat(); // ankur

        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(ThermoCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        ThermoCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        ThermoCloud(ThermoCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        ThermoCloud
        (
            const fvMesh& mesh,
            const word& name,
            const ThermoCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ThermoCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ThermoCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~ThermoCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const ThermoCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

            //- Return const access to thermo package
            inline const SLGThermo& thermo() const;

            //- Return const access to the carrier temperature field
            inline const volScalarField& T() const;

            //- Return const access to the carrier prressure field
            inline const volScalarField& p() const;


            // Sub-models

                //- Return reference to heat transfer model
                inline const HeatTransferModel<ThermoCloud<CloudType>>&
                    heatTransfer() const;


            // Integration schemes

                //-Return reference to velocity integration
                inline const scalarIntegrationScheme& TIntegrator() const;


            // Modelling options

                 //- Radiation flag
                inline bool radiation() const;

                 //- coupled radiation flag
                inline bool coupledRadiation() const; // ankur

                //- Return the number of bands
                inline label nBands() const; // ankur

                //- Radiation sum of parcel projected areas [m2]
                inline volScalarField::Internal& radAreaP(const label bandI); // ankur

                //- Radiation sum of parcel projected areas [m2]
                inline const volScalarField::Internal&
                    radAreaP(const label bandI) const; // ankur

                //- Radiation sum of parcel projected areas for scattering [m2]
                inline volScalarField::Internal& radAreaPSc(const label bandI); // ankur

                //- Radiation sum of parcel projected areas for scattering [m2]
                inline const volScalarField::Internal&
                    radAreaPSc(const label bandI) const; // ankur

                //- Radiation sum of parcel temperature^4 [K4]
                inline volScalarField::Internal& radT4(const label bandI); // ankur

                //- Radiation sum of parcel temperature^4 [K4]
                inline const volScalarField::Internal& radT4(const label bandI) const; // ankur

                //- Radiation sum of parcel projected area*temperature^4 [m2K4]
                inline volScalarField::Internal& radAreaPT4(const label bandI); // ankur

                //- Radiation sum of parcel temperature^4 [m2K4]
                inline const volScalarField::Internal&
                    radAreaPT4(const label bandI) const; // ankur

                //- Radiation sum of scattering * asyymetry factor
                inline volScalarField::Internal& radAreaPScAsy(const label bandI); // ankur

                //- Radiation sum of scattering * asyymetry factor
                inline const volScalarField::Internal& radAreaPScAsy(const label bandI) const; // ankur

                //- Radiation sum of normalization coefficients for in-scattering term
                inline volScalarField::Internal& CScat(const label bandI); // ankur

                //- Radiation sum of normalization coefficients for in-scattering term
                inline const volScalarField::Internal& CScat(const label bandI) const; // ankur

                //- Radiation dia-Banded model diameters
                inline const List<scalar>& diaVal() const; // ankur

            // Sources

                // Enthalpy

                    //- Sensible enthalpy transfer [J/kg]
                    inline volScalarField::Internal& hsTrans();

                    //- Sensible enthalpy transfer [J/kg]
                    inline const volScalarField::Internal&
                        hsTrans() const;

                    //- Return coefficient for carrier phase hs equation
                    inline volScalarField::Internal& hsCoeff();

                    //- Return const coefficient for carrier phase hs equation
                    inline const volScalarField::Internal&
                        hsCoeff() const;

                    //- Return sensible enthalpy source term [J/kg/m3/s]
                    inline tmp<fvScalarMatrix> Sh(volScalarField& hs) const;


                // Radiation - overrides thermoCloud virtual abstract members

                    //- Return tmp equivalent particulate emission
                    inline tmp<volScalarField> Ep(const label bandI) const; // ankur

                    //- Return tmp equivalent particulate absorption
                    inline tmp<volScalarField> ap(const label bandI) const; // ankur

                    //- Return tmp equivalent particulate scattering factor
                    inline tmp<volScalarField> sigmap(const label bandI) const; // ankur

                    //- Return tmp phase-function 
                    inline tmp<volScalarField> phaseFunc(const label bandI, const label sour, const label dest) const; // ankur

                 // Radiation - Access functions for radiation calculations for parcels

                     // - Update cloud level radiation properties in cells
                     inline void updateRadProperties(const label celli, const scalar ap, const label np, const scalar T4, const scalar dt); // ankur

                     // - Update cloud level radiation properties in cells, and compute quantities to calculate radiation source terms for parcels
                     inline List<scalar> radCalc(const label celli, const scalar ap, const scalar np, const scalar T4, const scalar dt); // ankur

                 // - Compute k*G for a parcel particle
                 inline scalar getkG(const scalar pDia) const;  // ankur

             // - Compute QEmm and emFrac of a parcle particle for all bands
             //inline scalar getEmm(const scalar pDia) const;
             inline void getQEmmFrac(const scalar pDia, const scalar pT, List<scalar>& QEmm, List<scalar>& emFrac) const; // ankur

             // - Compute Qabs of a parcle particle
             inline scalar getQAbs(const scalar pDia, const label bandI) const; // ankur

             // - Compute Qabs and Qsct of a parcle for all bands
             inline void getQAbsQSct(const scalar, const scalar nParticle, List<scalar>& QAbs, List<scalar>& QSct, List<scalar>& asyFac) const; // ankur

        // Check

            //- Maximum temperature
            inline scalar Tmax() const;

            //- Minimum temperature
            inline scalar Tmin() const;


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const ThermoCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve();

            //- Evolve the cloud
            void evolve();

            //- update properties of the cloud
            void updateProp(); // ankur


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

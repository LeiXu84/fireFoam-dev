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

Class
    Foam::pyroCUPOneDimV1

Description
    Reacting, 1-D pyrolysis model for CUP commodity

SourceFiles
    pyroCUPOneDimV1.C

\*---------------------------------------------------------------------------*/

#ifndef pyroCUPOneDimV1_H
#define pyroCUPOneDimV1_H

#include "reactingOneDim21CharOxi.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{


/*---------------------------------------------------------------------------*\
                      Class pyroCUPOneDimV1 Declaration
\*---------------------------------------------------------------------------*/

class pyroCUPOneDimV1
:
    public reactingOneDim21CharOxi
{
private:
 
            void initParams();    // Initializing the parameters for the CUP pyrolysis model

protected:

    // Protected data

            // Variables that could take different values at different boundary facets

            volScalarField tempUC_;      // temperature of the unit cell
            volScalarField massUC_;      // remaining mass of the unit cell
            volScalarField massUCCC_;      // remaining mass of CC in the unit cell
            volScalarField massUCPS_;      // remaining mass of PS in the unit cell
            volScalarField CpUC_;        // heat capacity of the unit cell
            volScalarField alphaUC_;     // absorptivity of the unit cell
            volScalarField emmUC_;       // emissivity of the unit cell
            volScalarField mdotUC_;      // mass loss rate from the outer unit cell
            volScalarField mdotIUC_;     // mass loss rate from inner unit cells  (from lumped model)
            volScalarField swccGone_;    // label to indicate if outer cardboard liner has burnt through
            volScalarField outerUCGone_; // label to indicate if outer unit cell has burnt through  
      
            volScalarField massInnerUC_;      // remaining representative mass of the inner unit cell on the boundary patch

            // Global variables for the CUP commodity (either each pallet load or all pallet loads).
            // Note that we will have to define global variables for each pallet load later on (for example,
            // mass remaining in each pallet load once the lumped model kicks in).

            dimensionedScalar heatGassUC_;      // Heat of gassification for outer unit cell
            dimensionedScalar heatGassUC2_;     // Second Heat of gassification for outer unit cell
            dimensionedScalar QFlameUC_;        // flame heat flux for outer unit cell
            dimensionedScalar QFlameUC2_;       // second regime flame heat flux for outer unit cell
            dimensionedScalar QEmmUC_;          // emissive heat flux from outer unit cell
            dimensionedScalar OLCrit_;          // OL criterion for it to exfoliate
            dimensionedScalar OUCEnSplit_;      //  Criterion for OUC to start sharing energy with IUC 
            dimensionedScalar multFacQFl_;      //  Multiplication factor for flame heat flux in case of concurrent burning of OUC and IUC or pure IUC burning

            dimensionedScalar massFracUC_;      // Mass-fraction of CC to be pyrolyzed in the first regime for outer unit cell
            bool TwoHeatGass_;                  // boolean for two heat of gassification approach (on or off)
            dimensionedScalar Tacti_;           // activation temperature for gassification reaction in the outer unit cell
            dimensionedScalar AGassReac_;       // Pre-exponential factor for the gassification reaction in the outer unit cell
            dimensionedScalar nGassReac_;       // reaction order for the gassification reaction in the outer unit cell
            dimensionedScalar initMassUC_;      // initial mass in the outer unit cell
            dimensionedScalar initMassUCCC_;    // initial mass of CC in the outer unit cell
            dimensionedScalar initMassUCPS_;    // initial mass of PS in the outer unit cell

            dimensionedScalar heatGassIUC_;   // Heat of gassification for the inner unit cells
            dimensionedScalar tempIUC_;        // Representative temperature for inner unit cells      
            dimensionedScalar emmIUC_;         // Representative emmissivity for inner unit cells      
            dimensionedScalar absIUC_;         // Representative absorptivity for inner unit cells      
            dimensionedScalar initMassIUC_;    // Initial mass for inner unit cells
            dimensionedScalar initMassIUCCC_;   // initial mass of CC in the inner unit cell
            dimensionedScalar initMassIUCPS_;   // initial mass of PS in the inner unit cell
            dimensionedScalar massIUC_;        // Remaining mass for inner unit cells (lumped total mass remaining for each pallet load)

            dimensionedScalar igniTempUC_;       // Ignition temperature for Constant-ignition-temperature pyrolysis model
            label           pModelTypeUC_;       // pyrolysis model for unit cell.. 1-> constant temp ignition, 2-> Arrhenius form

            dimensionedScalar hocPyrCC_;         // Representative heat of combustion for CC pyrolsate from outer unit cells
            dimensionedScalar hocPyrPS_;         // Representative heat of combustion for PS pyrolsate from outer unit cells
            dimensionedScalar hocPyrUC_;         // Representative heat of combustion for pyrolsate from outer unit cells
            dimensionedScalar hocPyrUC2_;        // Representative heat of combustion for pyrolsate from outer unit cells
            dimensionedScalar hocPyrIUC_;        // Representative heat of combustion for pyrolsate from inner unit cells

            // Variables to store convective and radiative heat flux values at the boundary patches
            volScalarField qradBnd_;      // radiative flux at the boundary
            volScalarField qconvBnd_;     // convective flux at the boundary
         
            // Vars to store the names of fuel_CC (corresponding to corrugated cardboard) and fuel_PS (polystyrene cup) species
            bool multiFuel_;             // boolean for multiFuel approach (on or off)
            word speciesCC_;             // fuel species for CC      
            word speciesPS_;             // fuel species for PS      


            // Variables that store fuel pallets configuration in a rack storage scenario

            label nPallets_;               // Number of fuel pallets

            label nXPallets_;              // Pallets in X-direction
            label nYPallets_;              // Pallets in Y-direction
            label nZPallets_;              // Pallets in Z-direction

            dimensionedScalar dXPallets_;  //  Distance between two pallets in X-direction (pallet center to pallet center)
            dimensionedScalar dYPallets_;  //  Distance between two pallets in Y-direction (center to center)
            dimensionedScalar dZPallets_;  //  Distance between two pallets in Z-direction (center to center)

            dimensionedScalar lXPallet_;  //  Pallet dimension in X-direction
            dimensionedScalar lYPallet_;  //  Pallet dimension in Y-direction
            dimensionedScalar lZPallet_;  //  Pallet dimension in Z-direction

            dimensionedVector botCorPalletOrigin_;  // Centroid for the pallet that is located at min X, min Y and min Z

            // Variables to track mass loss rates from individual pallets
            // This fields are of the size(# of pallets)

            scalarField massFluxOL_;    // mass loss rate from outer liner
            scalarField massFluxOUC_;   // mass loss rate from outer unit cell
            scalarField massFluxIUC_;   // mass loss rate from inner unit cell
            scalarField massFluxOCC_;    // mass loss rate for cardboard component of outer UC
            scalarField massFluxOPS_;    // mass loss rate for poly-styrene component of outer UC
            scalarField massFluxICC_;    // mass loss rate for cardboard component of inner UCs
            scalarField massFluxIPS_;    // mass loss rate for poly-styrene component of inner UCs

            // These arrays of size (# of pyrolysis patch faces)
            scalarField massFluxCCFracIUC_; // mass flux fraction of CC for inner unit cell
            scalarField massFluxCCFracOUC_; // mass flux fraction of CC for outer unit cell


            // Variables to store connectivity of patch faces to pallets

            labelList faceToPallet_;       //  Stores patch face to pallet # information
            labelList nPalletFaces_;       //  Stores number of patch faces on each pallet
            labelListList palletToFace_;  //  Stores pallet to patch face # information

            // Variables to track total mass (remaining and initial) in the pallets
            // This fields are of the size(# of pallets)

            scalarField totalMassIUC_;        // mass remaining in the inner unit cells for individual pallets
            scalarField diagTMIUC_;           // For diagnostics
            scalarField palletArea_;           // pallet area
            scalarField totalMassInitIUC_;    // Initial mass in the inner unit cells for individual pallets
            scalarField totalMassOUC_;        // mass remaining in the outer unit cells for individual pallets
            scalarField diagTMOUC_;           // For diagnostics
            scalarField totalMassInitOUC_;    // Initial mass in the outer unit cells for individual pallets
            scalarField totalMassOL_;         // mass remaining in the outer liner for individual pallets
            scalarField diagTMOL_;            // For diagnostics
            scalarField totalMassInitOL_;     // Initial mass in the outer liner for individual pallets

            // Variables for back boundary of OL - OUC coupling

            label     numCoupledPatches_;      // number of couples patches
            labelList backBndPatch_;          // Size (# of coupled patches)...  Stores back bnd patch number

            // Variables and mapping arrays defined for improving and streamlining code computations..

            labelList coupledPatchID_;        // Array to store patch numbers
            labelList startFaceID_;           // Starting face index for a patch in the consolidated array for coupled faces
            labelListList numPatchFaceCt_;    // Number of coupled patch faces corresponding to OL, OUC, OUC2, and IUC regimes             
  
            labelListList facesOL_;           // Mapping array for OL faces 
            labelListList facesOUC_;          // Mapping array for OUC faces  - regime 1
            labelListList facesOUC2_;         // Mapping array for OUC faces  - regime 2
            labelListList facesIUC_;          // Mapping array for IUC faces

            scalarField pyrHOC_;               // scalar field to store heat of combustion for pyrolysis faces 
            labelList   patchIDToCoupledID;   //  Mapping array to get coupled patch ID for a patch.. (will return valid coupled patch ID for a coupled patch otherwise return -1 for non-coupled patches) 
            //labelList   firstCallSetHOC_;      // Index to store whether this is the first call to getPyroHOC function or not (1-> first call, 0-> not a first call)
                                               // This list is of size (# of couples patches)

            // Indexes (that represent various regimes) to allow for operating on entire fields (instead of facet by facet on a coupled patch) 

            volScalarField IOL_;         // Index for the outer liner
            volScalarField IOUCHU_;      // Index for the outer unit cell- heat up period
            volScalarField IOUCR1_;      // Index for the outer unit cell- regime 1
            volScalarField IOUCR2_;      // Index for the outer unit cell- regime 2
            volScalarField IIUC_;        // Index for the inner unit cell


            // Initial mass and aggregated mass loss from outer liner (to allow for determining if the outer liner is burnt by comparing aggreated mass loss with the initial mass)
                      
            volScalarField mOLInit_;         // Initial mass in the OL liner (this is not the mass density but the total virgin mass in that 1-D section)
            volScalarField mOLLost_;         // Aggregated lost mass in the OL liner (not the lost mass density but the total lost virgin mass from that 1-D section)

            // Variables defined to store emissivity and absorptivity values for different stages of outer unit cell
            dimensionedScalar emmIOUCHU_;         // Representative emmissivity for outer unit cell- heat up period     
            dimensionedScalar absIOUCHU_;         // Representative absorptivity for outer unit cell- heat up period
            dimensionedScalar emmIOUCR1_;         // Representative emmissivity for outer unit cell- regime 1..    
            dimensionedScalar absIOUCR1_;         // Representative absorptivity for outer unit cell- regime 1..
            dimensionedScalar emmIOUCR2_;         // Representative emmissivity for outer unit cell- regime 2..     
            dimensionedScalar absIOUCR2_;         // Representative absorptivity for outer unit cell- regime 2..

            // variables defined to compute convective and enclosure losses as fractions of incident radiation 
            dimensionedScalar UCEnLossFr_;        // Enclosure loss fraction
            bool EnLossFracSpec_;                 // boolean for enclosure  loss specification
            dimensionedScalar UCEnLossFixed_;     // Enclosure loss fixed value

            dimensionedScalar UCConvLossFr_;      // Convective loss fraction
            dimensionedScalar UCConvLossFrR2_;    // Convective loss fraction
            bool ConvLossFracSpec_;               // boolean for convective  loss specification
            dimensionedScalar UCConvLossFixed_;   // convective loss fixed value

            // variable defined to store the extra contribution from flame in large fire scenarios 
            dimensionedScalar QFlameExtra_;       // second regime flame heat flux for outer unit cell
            dimensionedScalar QFlameExtraOUCR1_;  // second regime flame heat flux for outer unit cell

            // variable to store energy releases inside UC due to burning of fuel
            volScalarField enerRelUC_;            // Energy releases inside the unit cell

            // Variables to store fields for heat of combustion and face to pallet ordering
         
            volScalarField pyrolHOC_;             // Heat of combustion for the pyrolysate
            PtrList<volScalarField> face2Pallet_;   // volScalarFields for faceToPallet mapping  (one volcalarField for each pallet- stored 0 or 1 at patch boundaries indicating whether face belongs to this pallet or not) .. This is to enable field-wise operations..

    // Protected member functions

        //- Update submodels
        virtual void updateFields();

        //- Update emissivity and absorptivity at boundaries
        virtual void updateBndEmmAbs();

        //- Solve for species equations
        virtual void solveSpeciesMass();

        //- Update the BCs for fuel species if multiFuel option active
        virtual void updateFuelFluxMassFrac();

        // Equations

            //- Solve energy
            //virtual void solveEnergy();

public:

      // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the pyrolysis equations
            virtual void evolveRegion();

            //- Return the region absorptivity [1/m]
            virtual tmp<volScalarField> kappaRad() const;

            //- Public member functions to set the qradBnd and qconvBnd fieldss
            void setqrad(const scalarField& qradin, label patchI);
            void setQconv(const scalarField& Qconvin, label patchI);

            //- Public member function to get the HOC for pyrolsate from the pyrolysis model
            void getPyroHOC(scalarField& pyrHoc, const scalar& pyrHocCC, label patchI);
            

    //- Runtime type information
    TypeName("pyroCUPOneDimV1");


    // Constructors

        //- Construct from type name and mesh
        pyroCUPOneDimV1
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name, mesh and dictionary
        pyroCUPOneDimV1
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );


    //- Destructor
    virtual ~pyroCUPOneDimV1();

    // I-O

    //- Provide some feedback
    virtual void info() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "pyroCUPOneDimV1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

#include "pyroCUPOneDimV1.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvMatrices.H"
#include "absorptionEmissionModel.H"
#include "fvcLaplacian.H"
#include "physicoChemicalConstants.H"
#include "nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField.H"
#include "constHTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pyroCUPOneDimV1, 0);

addToRunTimeSelectionTable(pyrolysisModel, pyroCUPOneDimV1, mesh);
addToRunTimeSelectionTable(pyrolysisModel, pyroCUPOneDimV1, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void pyroCUPOneDimV1::initParams()
{
      igniTempUC_.value()  = coeffs().lookupOrDefault<scalar>("igniTempUC",800.);
      scalar CpUCVal       = coeffs().lookupOrDefault<scalar>("CpUC",600.);
      scalar emmUCVal      = coeffs().lookupOrDefault<scalar>("emissivityUC",0.6);
      scalar absUCVal      = coeffs().lookupOrDefault<scalar>("absorptivityUC",0.75);


      emmIOUCHU_.value() = emmUCVal; 
      absIOUCHU_.value() = absUCVal; 
      emmIOUCR1_.value()   = coeffs().lookupOrDefault<scalar>("emissivityUCR1",emmUCVal);
      emmIOUCR2_.value()   = coeffs().lookupOrDefault<scalar>("emissivityUCR2",emmUCVal);
      absIOUCR1_.value()   = coeffs().lookupOrDefault<scalar>("absorptivityUCR1",absUCVal);
      absIOUCR2_.value()   = coeffs().lookupOrDefault<scalar>("absorptivityUCR2",absUCVal);
      //emmIOUCR1_.value() = emmUCVal; 
      //emmIOUCR2_.value() = emmUCVal; 
      //absIOUCR1_.value() = absUCVal; 
      //absIOUCR2_.value() = absUCVal; 

      heatGassUC_.value()  = coeffs().lookupOrDefault<scalar>("heatGassUC",1.8e+6);
      heatGassUC2_.value() = coeffs().lookupOrDefault<scalar>("heatGassUC2",6e+6);


      QFlameUC_.value()  = coeffs().lookupOrDefault<scalar>("QFlameUC",30000.);
      QFlameUC2_.value()  = coeffs().lookupOrDefault<scalar>("QFlameUC2",36000.);
      QFlameExtra_.value()  = coeffs().lookupOrDefault<scalar>("QFlameExtra",10000.);
      QFlameExtraOUCR1_.value()  = coeffs().lookupOrDefault<scalar>("QFlameExtraOUCR1",10000.);
      QEmmUC_.value()  = coeffs().lookupOrDefault<scalar>("QEmmUC",15000.);
      OLCrit_.value()  = coeffs().lookupOrDefault<scalar>("OLCrit",0.04);
      OUCEnSplit_.value()  = coeffs().lookupOrDefault<scalar>("OUCEnSplit",0.02);
      multFacQFl_.value()  = coeffs().lookupOrDefault<scalar>("multFacQFl",1.);

      EnLossFracSpec_        = coeffs().lookupOrDefault<bool>("EnLossFracSpec",true);  
      UCEnLossFr_.value()    = coeffs().lookupOrDefault<scalar>("UCEnLossFr",0.2);
      UCEnLossFixed_.value() = coeffs().lookupOrDefault<scalar>("UCEnLossFixed",20000.);

      ConvLossFracSpec_        = coeffs().lookupOrDefault<bool>("ConvLossFracSpec",true);  
      UCConvLossFr_.value()    = coeffs().lookupOrDefault<scalar>("UCConvLossFr",0.1);
      UCConvLossFrR2_.value()  = coeffs().lookupOrDefault<scalar>("UCConvLossFrR2",0.1);
      UCConvLossFixed_.value() = coeffs().lookupOrDefault<scalar>("UCConvLossFixed",10000.);

      massFracUC_.value()  = coeffs().lookupOrDefault<scalar>("CCMassFracUC",0.6);
      TwoHeatGass_         = coeffs().lookupOrDefault<bool>("TwoHeatGass",false);  
      heatGassIUC_.value() = coeffs().lookupOrDefault<scalar>("heatGassInnerUC",3e+6);

      initMassUCCC_.value()  = coeffs().lookupOrDefault<scalar>("initMassUCCC",2.2);
      initMassUCPS_.value()  = coeffs().lookupOrDefault<scalar>("initMassUCPS",3.2);

      initMassIUCCC_.value()  = coeffs().lookupOrDefault<scalar>("initMassIUCCC",2.2);
      initMassIUCPS_.value()  = coeffs().lookupOrDefault<scalar>("initMassIUCPS",3.2);

      // In the pyro model, the initial mass for UC and IUC are computed from the initial masses for CC and PS
      initMassUC_.value()  = initMassUCCC_.value() + initMassUCPS_.value();  
      initMassIUC_.value() = initMassIUCCC_.value() + initMassIUCPS_.value();  

      hocPyrCC_.value()    = coeffs().lookupOrDefault<scalar>("hocPyrCC",1.3e+7);    
      hocPyrPS_.value()    = coeffs().lookupOrDefault<scalar>("hocPyrPS",2.6e+7);    
      // In the pyro model, the hoc is computed from hoc for CC and PS, and the inital masses..


      hocPyrUC_.value() = hocPyrCC_.value(); 

      hocPyrUC2_.value() = 
               ((1. - massFracUC_.value())*initMassUCCC_.value()*hocPyrCC_.value() + initMassUCPS_.value()*hocPyrPS_.value()) /
               ((1. - massFracUC_.value())*initMassUCCC_.value() + initMassUCPS_.value());

      hocPyrIUC_.value() = 
               (initMassIUCCC_.value()*hocPyrCC_.value() + initMassIUCPS_.value()*hocPyrPS_.value()) /
               (initMassIUCCC_.value() + initMassIUCPS_.value());


      tempIUC_.value()     = coeffs().lookupOrDefault<scalar>("TempInnerUC",800.);    
      emmIUC_.value()      = coeffs().lookupOrDefault<scalar>("emissivityInnerUC",0.6);    
      absIUC_.value()      = coeffs().lookupOrDefault<scalar>("absorptivityInnerUC",0.6);    
     
      multiFuel_           = coeffs().lookupOrDefault<bool>("multiFuel",false);  
   
      if (multiFuel_)
      {
         speciesCC_           = coeffs().lookupOrDefault<word>("speciesCC","none");  
         speciesPS_           = coeffs().lookupOrDefault<word>("speciesPS","none");  

        if (speciesCC_ == "none" || speciesPS_ == "none")
        {
           FatalErrorIn
           (
            "pyroCUPOneDimV1::initParams"
            "("
            ")"
           )
            << " SpeciesCC and SpeciesPS need to be specified in the dictionary.\n"
            << nl << nl 
            << coeffs() << nl << nl
            << " Please update the dictionary with the required input." << exit(FatalError);             
        }
      }

      Info << "Igni Temp: " << igniTempUC_.value() << endl; 
      Info << "Cp UC:     " << CpUCVal << endl;
      Info << "emm UC:    " << emmUCVal << endl;
      Info << "abs UC:    " << absUCVal << endl;

      Info << "OL Crit:    " << OLCrit_.value() << endl;
      Info << "OUC Energy Split Crit:    " << OUCEnSplit_.value() << endl;
      Info << "Flame heat flux multiplication factor during concurrent burning of OUC and IUC or during pure IUC burning:  " << multFacQFl_.value() << endl;

      if (EnLossFracSpec_)
      {
         Info << "Enclosure loss fraction:    " << UCEnLossFr_.value() << endl;
      }
      else
      {
         Info << "Enclosure loss Flux:    " << UCEnLossFixed_.value() << endl;
      }

      if (ConvLossFracSpec_)
      {
         Info << "Convecive loss fraction:    " << UCConvLossFr_.value() << endl;
      }
      else
      {
         Info << "Convective loss flux:   " << UCConvLossFixed_.value() << endl;
      }

      Info << "emm UCHU:    " << emmIOUCHU_.value() << endl;
      Info << "abs UCHU:    " << absIOUCHU_.value() << endl;
      Info << "emm UCR1:    " << emmIOUCR1_.value() << endl;
      Info << "abs UCR1:    " << absIOUCR1_.value() << endl;
      Info << "emm UCR2:    " << emmIOUCR2_.value() << endl;
      Info << "abs UCR2:    " << absIOUCR2_.value() << endl;

      Info <<"heatGass UC:" << heatGassUC_.value() << endl;
      if (TwoHeatGass_)
      {
          Info <<"heatGass2 UC:" << heatGassUC2_.value() << endl;
      }
      Info <<"heatGassIUC:" << heatGassIUC_.value() << endl;
      Info <<"Mass/A UC:  " << initMassUC_.value() << endl;
      Info <<"Mass/A UC CC:  " << initMassUCCC_.value() << endl;
      Info <<"Mass/A UC PS:  " << initMassUCPS_.value() << endl;
      Info <<"Mass/A IUC: " << initMassIUC_.value() << endl;
      Info <<"hocPyr UC:  " << hocPyrUC_.value() << endl;
      Info <<"hocPyr IUC: " << hocPyrIUC_.value() << endl;
      Info <<"Temp IUC:   " << tempIUC_.value() << endl;
      Info <<"emm IUC:    " << emmIUC_.value() << endl;
      Info <<"abs IUC:    " << absIUC_.value() << endl;
      Info <<"multiFuel:  " << multiFuel_ << endl;
      if (multiFuel_)
      {
        Info <<"speciesCC:  " << speciesCC_ << endl;
        Info <<"speciesPS:  " << speciesPS_ << endl;
      }


        forAll(intCoupledPatchIDs_, i)
        {
          const label patchI = intCoupledPatchIDs_[i];
          scalarField& Cp_UC    = CpUC_.boundaryFieldRef()[patchI];
          scalarField& alpha_UC = alphaUC_.boundaryFieldRef()[patchI];
          scalarField& emm_UC   = emmUC_.boundaryFieldRef()[patchI];
   
          emm_UC   = emmUCVal; // 0.6;
          alpha_UC = absUCVal; // 0.75; // 0.65;

          Cp_UC    = CpUCVal; // 600.; // 1000.;     
        }


      // Setting variables for fuel pallets configuration

      nXPallets_  = coeffs().lookupOrDefault<label>("nXPallets",2);
      nYPallets_  = coeffs().lookupOrDefault<label>("nYPallets",2);
      nZPallets_  = coeffs().lookupOrDefault<label>("nZPallets",2);

      nPallets_ = nXPallets_ * nYPallets_ * nZPallets_;

      dXPallets_.value()  = coeffs().lookupOrDefault<scalar>("dXPallets",1.2192);
      dYPallets_.value()  = coeffs().lookupOrDefault<scalar>("dYPallets",1.2192);
      dZPallets_.value()  = coeffs().lookupOrDefault<scalar>("dZPallets",1.524);
 
      lXPallet_.value()   = coeffs().lookupOrDefault<scalar>("lXPallet",1.0668);
      lYPallet_.value()   = coeffs().lookupOrDefault<scalar>("lYPallet",1.0668);
      lZPallet_.value()   = coeffs().lookupOrDefault<scalar>("lZPallet",1.0668);

      Info << " dX dY dZ: " << dXPallets_.value() << " " << dYPallets_.value() << " " << dZPallets_.value() << endl;
      Info << " lX lY lZ: " << lXPallet_.value() << " " << lYPallet_.value() << " " << lZPallet_.value() << endl;

      botCorPalletOrigin_.value() = coeffs().lookupOrDefault<vector>("botCorPalletOrigin", Foam::vector(-0.6096,-0.6096,1.143));

      Info << "botCor Origin: " << botCorPalletOrigin_.value() << endl;

      // Settings variables that track mass loss

      massFluxOL_      = scalarField(nPallets_,0.);
      massFluxOUC_     = scalarField(nPallets_,0.);
      massFluxIUC_     = scalarField(nPallets_,0.);
      massFluxOCC_     = scalarField(nPallets_,0.);
      massFluxOPS_     = scalarField(nPallets_,0.);
      massFluxICC_     = scalarField(nPallets_,0.);
      massFluxIPS_     = scalarField(nPallets_,0.);


      // Settings variables that store connectivity between patch faces and pallet #s

      label pyrolysisFacesCount = 0;

      numCoupledPatches_ = 0;

      forAll(intCoupledPatchIDs_, i)
      {
          const label patchI = intCoupledPatchIDs_[i];
          pyrolysisFacesCount = pyrolysisFacesCount + regionMesh().boundaryMesh()[patchI].size();

          numCoupledPatches_ = numCoupledPatches_ + 1;   // Used in OL-OUC coupling
      }
      backBndPatch_ = labelList(numCoupledPatches_,-1);  // Used in OL_OUC coupling 

      massFluxCCFracIUC_ = scalarField(pyrolysisFacesCount,0.);
      massFluxCCFracOUC_ = scalarField(pyrolysisFacesCount,0.);


      faceToPallet_    = labelList(pyrolysisFacesCount,-1);
      Info << "Pyro face Ct: " << pyrolysisFacesCount << endl;

      nPalletFaces_    = labelList(nPallets_,0);

      point minPt(0,0,0);
      point maxPt(1,1,1);
      List<boundBox> boundBoxPallets(nPallets_,boundBox(minPt,maxPt));


      forAll(boundBoxPallets,palletID)
      {
          point& minPtPallet=boundBoxPallets[palletID].min();
          point& maxPtPallet=boundBoxPallets[palletID].max();

          label iX,iY,iZ,remY,remZ,planarID;

          iZ = 1 + (palletID / (nXPallets_ * nYPallets_));
          remZ = (palletID+1) - (iZ-1) * nXPallets_ * nYPallets_;

          if (remZ == 0)
          {
             planarID =  nXPallets_ * nYPallets_;
          }
          else
          {
             planarID = remZ;
          }

          iY   = 1 + ((planarID-1) / nXPallets_);
          remY = planarID - (iY-1) * nXPallets_;

          if (remY == 0)
          {
             iX =  nXPallets_;
          }
          else
          {
             iX = remY;
          }


          minPtPallet.x()= botCorPalletOrigin_.value().x() + (iX-1)*dXPallets_.value() - (0.5+0.02)*lXPallet_.value();
          maxPtPallet.x()= botCorPalletOrigin_.value().x() + (iX-1)*dXPallets_.value() + (0.5+0.02)*lXPallet_.value();

          minPtPallet.y()= botCorPalletOrigin_.value().y() + (iY-1)*dYPallets_.value() - (0.5+0.02)*lYPallet_.value();
          maxPtPallet.y()= botCorPalletOrigin_.value().y() + (iY-1)*dYPallets_.value() + (0.5+0.02)*lYPallet_.value();

          minPtPallet.z()= botCorPalletOrigin_.value().z() + (iZ-1)*dZPallets_.value() - (0.5+0.02)*lZPallet_.value();
          maxPtPallet.z()= botCorPalletOrigin_.value().z() + (iZ-1)*dZPallets_.value() + (0.5+0.02)*lZPallet_.value();

         Info << "pallet " << palletID << " ix: " << iX << " iy: " << iY << " iz: " << iZ << nl;
         Info << "pallet box " << palletID << " " << boundBoxPallets[palletID].min() << " " << boundBoxPallets[palletID].max() << nl;
      }

// Setting volScalarFields for face to pallet mapping

        face2Pallet_.setSize(nPallets_);

        forAll(nPalletFaces_, i)
        {
            face2Pallet_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "f2Pallet_" + Foam::name(i),
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar("zero", dimless, 0.)
                )
            );
        }
// 

      forAll(intCoupledPatchIDs_, i)
      {

          const label patchI = intCoupledPatchIDs_[i];          
          const polyPatch&  pPatch = regionMesh().boundaryMesh()[patchI];

          label pyroFaceID = 0;

          label pID = -1;
          forAll(pPatch,faceI)
          {
            const point& fCenter = pPatch.faceCentres()[faceI];

            forAll(boundBoxPallets,palletID)
            {
                // checking if point (face center lies inside the bounding box for this pallet....
                if (boundBoxPallets[palletID].contains(fCenter))
                {
                    faceToPallet_[pyroFaceID]=palletID;
                    nPalletFaces_[palletID]=nPalletFaces_[palletID]+1;
                    break;
                }
            }
       
             if (faceToPallet_[pyroFaceID] == -1)
             {
                 Info << "Untracked: " << fCenter << endl; 
             }

             pID = faceToPallet_[pyroFaceID];
             if (pID > -1) face2Pallet_[pID].boundaryFieldRef()[patchI][faceI] = 1.;   // setting this to 1 for the pallet this face belongs to..  

             pyroFaceID = pyroFaceID + 1;
          }

      }


      forAll(faceToPallet_,faceID)
      {
          if (faceToPallet_[faceID] == -1)
          {
              FatalErrorIn
              (
               "pyroCUPOneDimV1::initParams"
               "("
               ")"
              )
               << " Error in setting the pyrolysis patch faces to pallet ordering.\n"
               << nl << nl
               << " Please check the inputs in pyrolysis dictionary."
               << nl << nl
               << coeffs() << nl << nl
               << exit(FatalError);
          }
      }


      palletToFace_    =  labelListList(nPallets_);

      forAll(palletToFace_,lListID)
      {
          palletToFace_[lListID].setSize(nPalletFaces_[lListID]);
      }

      labelList palletsCt(nPallets_,0);
      forAll(faceToPallet_,faceID)
      {
          label palletID = faceToPallet_[faceID];
          palletToFace_[palletID][palletsCt[palletID]] = faceID;
          palletsCt[palletID] = palletsCt[palletID] + 1;
      }


      // Settings variables that track the remaining mass in outer liners, outer UC and inner unit cells
      // for each pallets..
      // Initial mass of the OL is also set in this routine...

      totalMassIUC_       =   scalarField(nPallets_,0.);
      palletArea_         =   scalarField(nPallets_,0.);
      totalMassInitIUC_   =   scalarField(nPallets_,0.);
      totalMassOUC_       =   scalarField(nPallets_,0.);
      totalMassInitOUC_   =   scalarField(nPallets_,0.);
      totalMassOL_        =   scalarField(nPallets_,0.);
      totalMassInitOL_    =   scalarField(nPallets_,0.);

      diagTMIUC_          =   scalarField(nPallets_,0.);
      diagTMOUC_          =   scalarField(nPallets_,0.);
      diagTMOL_           =   scalarField(nPallets_,0.);


      const label vIndex = solidThermo_.composition().species()["v"];
      const label chIndex = solidThermo_.composition().species()["char"];

      scalar rhoV = solidThermo_.composition().rho(vIndex,  10000, 298);
      scalar rhoC = solidThermo_.composition().rho(chIndex, 10000, 298);
      scalar vToGasFac = (rhoV-rhoC)/rhoV; //virgin mass to vapor pyrolsate conversion (i.e., 1 gm of virgin leads to this amount of vapor pyrolsate)..


      label localPyrolysisFaceI = 0;
      forAll(intCoupledPatchIDs_, i)
      {
          const label patchI = intCoupledPatchIDs_[i];
          const scalarField& mUC = massUC_.boundaryField()[patchI];
          const scalarField& mIUC = massInnerUC_.boundaryField()[patchI];

          scalarField& mOLInit_patch   = mOLInit_.boundaryFieldRef()[patchI];
          scalarField& mOLLost_patch   = mOLLost_.boundaryFieldRef()[patchI];

          const scalarField& cellV = regionMesh().V();

          const scalarField& fArea  = regionMesh().boundary()[patchI].magSf();

          forAll(mUC, faceI)
          {
              label palletID = faceToPallet_[localPyrolysisFaceI];

              const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI];
              scalar vMass     = 0.0;
              scalar vIniMass = 0.0;
              forAll(cells, k)
              {
                  const label cellI = cells[k];
                  vMass     += Ys_[vIndex][cellI] * rho_[cellI] *cellV[cellI];
                  vIniMass += rhoV * cellV[cellI];
              }
      
              mOLInit_patch[faceI] = vIniMass;          // Assigning the virgin mass (initial, at t=0) to the mOLInit_ field  
              mOLLost_patch[faceI] = (vIniMass - vMass)*vToGasFac;  // Assigning the pyrolsate generated from lost virgin mass to the mOLLost_ field  

              palletArea_[palletID]      = palletArea_[palletID] + fArea[faceI];

              totalMassOL_[palletID]     = totalMassOL_[palletID] + vMass;
              totalMassInitOL_[palletID] = totalMassInitOL_[palletID] + vIniMass;

              totalMassOUC_[palletID]      = totalMassOUC_[palletID] + mUC[faceI]*fArea[faceI];
              totalMassInitOUC_[palletID]  = totalMassInitOUC_[palletID] + initMassUC_.value()*fArea[faceI];

              totalMassIUC_[palletID]      = totalMassIUC_[palletID] + mIUC[faceI]*fArea[faceI];
              totalMassInitIUC_[palletID]  = totalMassInitIUC_[palletID] + initMassIUC_.value()*fArea[faceI];

              localPyrolysisFaceI++;
          }
      }


// ----
// Identifying patches on the back boundary for the outer liner; for the OL-OUC coupling 

      const fvPatchList& meshPatches = regionMesh().boundary();

      volScalarField& T1 = solidThermo_.T();
      const volScalarField::Boundary& Tbf = T1.boundaryField();

      localPyrolysisFaceI = 0;

      forAll(intCoupledPatchIDs_, i)
      {
          const label patchI = intCoupledPatchIDs_[i];

          forAll(meshPatches,patchj)
          {
            const fvPatchScalarField& Tp = Tbf[patchj];

            if (isA<Foam::constHTemperatureFvPatchScalarField>(Tp))
            {
              if ( Tp.patch().start() ==  boundaryFaceOppositeFace_[localPyrolysisFaceI])
              {
               backBndPatch_[i] = patchj; 
              } 
            }

          } 

          localPyrolysisFaceI = localPyrolysisFaceI + regionMesh().boundaryMesh()[patchI].size();

      }

// ------------
     // Variables and mapping arrays for improving code performance
     //  Indexes for allowing field-wise operations are also initialized here at the start/restart of the run..

     coupledPatchID_ = labelList(numCoupledPatches_,-1);
     startFaceID_    = labelList(numCoupledPatches_,-1);
        Info << "startFaceID_: " << startFaceID_ << endl;

     numPatchFaceCt_  =  labelListList(numCoupledPatches_);
     facesOL_         =  labelListList(numCoupledPatches_);
     facesOUC_        =  labelListList(numCoupledPatches_);
     facesOUC2_       =  labelListList(numCoupledPatches_);
     facesIUC_        =  labelListList(numCoupledPatches_);

      patchIDToCoupledID = labelList(regionMesh().boundaryMesh().types().size(),-1);
     pyrolysisFacesCount = 0; 
     forAll(intCoupledPatchIDs_, i)
     {
          const label patchI = intCoupledPatchIDs_[i];
          pyrolysisFacesCount = pyrolysisFacesCount + regionMesh().boundaryMesh()[patchI].size();
       
          patchIDToCoupledID[patchI] = i;
     }
     pyrHOC_ = scalarField(pyrolysisFacesCount,0.);
     
 
     // Mass threshold to decide on transition from regime 1 to regime 2
     scalar thresholdUCMass = initMassUC_.value()-massFracUC_.value()*initMassUCCC_.value();

     pyrolysisFacesCount = 0; 

     localPyrolysisFaceI = 0;

     forAll(intCoupledPatchIDs_, i)
     {
        const label patchI = intCoupledPatchIDs_[i];
        coupledPatchID_[i] = patchI; 
      
        startFaceID_[i] = pyrolysisFacesCount;

        pyrolysisFacesCount = pyrolysisFacesCount + regionMesh().boundaryMesh()[patchI].size();

        numPatchFaceCt_[i].setSize(4);
        numPatchFaceCt_[i] = 0;
       
        label faceCt=regionMesh().boundaryMesh()[patchI].size();
        facesOL_[i].setSize(faceCt);
        facesOUC_[i].setSize(faceCt);
        facesOUC2_[i].setSize(faceCt);
        facesIUC_[i].setSize(faceCt);

        facesOL_[i] = -1;
        facesOUC_[i] = -1;
        facesOUC2_[i] = -1;
        facesIUC_[i] = -1;

        const scalarField& ccBurntSt = swccGone_.boundaryField()[patchI];
        const scalarField& UCBurntSt = outerUCGone_.boundaryField()[patchI];
        const scalarField& m_UC      = massUC_.boundaryField()[patchI];           


        scalarField& IOL_patch       = IOL_.boundaryFieldRef()[patchI];
        scalarField& IOUCHU_patch    = IOUCHU_.boundaryFieldRef()[patchI];
        scalarField& IOUCR1_patch    = IOUCR1_.boundaryFieldRef()[patchI];
        scalarField& IOUCR2_patch    = IOUCR2_.boundaryFieldRef()[patchI];
        scalarField& IIUC_patch      = IIUC_.boundaryFieldRef()[patchI];

        scalarField& pyrHOC_patch    = pyrolHOC_.boundaryFieldRef()[patchI];

        const scalarField& T_UC      = tempUC_.boundaryField()[patchI];           

        forAll(ccBurntSt,faceI)
        {
            if (ccBurntSt[faceI] > 0.)  // unburnt val is -1.0 and burnt val is 1.0
            {
               if (UCBurntSt[faceI] > 0.)
               {
                   facesIUC_[i][numPatchFaceCt_[i][3]] = faceI;
                   numPatchFaceCt_[i][3] = numPatchFaceCt_[i][3] + 1;
                   pyrHOC_[localPyrolysisFaceI] = hocPyrIUC_.value();

                   IIUC_patch[faceI] = 1.;
                   pyrHOC_patch[faceI] = hocPyrIUC_.value();
               }
               else
               {
                 if (T_UC[faceI] < igniTempUC_.value()) 
                 {

                   IOUCHU_patch[faceI] = 1.;
                 }
                 else 
                 { 
                   // Things to do when outer unit cell is burning
                   if (TwoHeatGass_)
                   {  
                       if ( m_UC[faceI] > thresholdUCMass) 
                       {
                           facesOUC_[i][numPatchFaceCt_[i][1]] = faceI;
                           numPatchFaceCt_[i][1] = numPatchFaceCt_[i][1] + 1;
                           pyrHOC_[localPyrolysisFaceI] = hocPyrUC_.value();

                          IOUCR1_patch[faceI] = 1.;
                          pyrHOC_patch[faceI] = hocPyrUC_.value();
                       }
                       else
                       {
                           facesOUC2_[i][numPatchFaceCt_[i][2]] = faceI;
                           numPatchFaceCt_[i][2] = numPatchFaceCt_[i][2] + 1;
                           pyrHOC_[localPyrolysisFaceI] = hocPyrUC2_.value();

                          IOUCR2_patch[faceI] = 1.;
                          pyrHOC_patch[faceI] = hocPyrUC2_.value();
                       }

                   }
                   else
                   {
                       facesOUC_[i][numPatchFaceCt_[i][1]] = faceI;
                       numPatchFaceCt_[i][1] = numPatchFaceCt_[i][1] + 1;
                       pyrHOC_[localPyrolysisFaceI] = hocPyrUC_.value();

                       IOUCR1_patch[faceI] = 1.;
                       pyrHOC_patch[faceI] = hocPyrUC_.value();
                   }
                 } 

               }
           }
           else
           {
               facesOL_[i][numPatchFaceCt_[i][0]] = faceI;
               numPatchFaceCt_[i][0] = numPatchFaceCt_[i][0] + 1;

               IOL_patch[faceI] = 1.;
           }
           localPyrolysisFaceI++;
       }


     }


        Info << "startFaceID_: " << startFaceID_ << endl;
        Info << "coupledPatchID_: " << coupledPatchID_ << endl;

//---------

}

void pyroCUPOneDimV1::updateFields()
{
    reactingOneDim21CharOxi::updateFields();

}



void pyroCUPOneDimV1::updateBndEmmAbs()
{
/*
// Content of this function has been removed since this function is not used anymore.
*/
}


void pyroCUPOneDimV1::getPyroHOC(scalarField& pyrHoc, const scalar& hocPyrCC, label patchI)
{
        const scalarField& ccBurntSt      = swccGone_.boundaryField()[patchI];
        //const scalarField& UCBurntSt      = outerUCGone_.boundaryField()[patchI];
        //const scalarField& m_UC           = massUC_.boundaryField()[patchI];           
        const scalarField& pyrHOC_patch   = pyrolHOC_.boundaryField()[patchI];
  
        // Implementing new algorithm for improving code performance
        //label coupledPatchIndex = patchIDToCoupledID[patchI]; 
        //label startFaceIndex    = startFaceID_[coupledPatchIndex];        

        //label patchSize = m_UC.size();
        pyrHoc  = 0.5*(1.+ccBurntSt)*pyrHOC_patch + 0.5*(1.-ccBurntSt)*hocPyrCC;

}


void pyroCUPOneDimV1::updateFuelFluxMassFrac()
{

    const fvMesh& gasMesh = this->primaryMesh();
    const volScalarField& fuelCC = gasMesh.lookupObject<volScalarField>(speciesCC_); 
    const volScalarField::Boundary& fuelCC_bf = fuelCC.boundaryField();

    const volScalarField& fuelPS = gasMesh.lookupObject<volScalarField>(speciesPS_); 
    const volScalarField::Boundary& fuelPS_bf = fuelPS.boundaryField();

    massFluxCCFracIUC_ = 1.;
    massFluxCCFracOUC_ = 1.;

    label localPyrolysisFaceI = 0;

    // Mass threshold to decide on transition from regime 1 to regime 2
    scalar thresholdUCMass = initMassUC_.value()-massFracUC_.value()*initMassUCCC_.value();
    
    scalar UC_CCFrac = (1.-massFracUC_.value())*initMassUCCC_.value()*hocPyrCC_.value() /
                   ( (1.-massFracUC_.value())*initMassUCCC_.value()*hocPyrCC_.value() + initMassUCPS_.value()*hocPyrPS_.value() );

    scalar IUC_CCFrac = initMassIUCCC_.value()*hocPyrCC_.value() /
                   ( initMassIUCCC_.value()*hocPyrCC_.value() + initMassIUCPS_.value()*hocPyrPS_.value() );


    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        const scalarField& ccBurntSt = swccGone_.boundaryField()[patchI];
        const scalarField& UCBurntSt = outerUCGone_.boundaryField()[patchI];
        scalarField& m_UC     = massUC_.boundaryFieldRef()[patchI];           

        scalarField massFluxFrac_CC(ccBurntSt.size(),1.); 
        scalarField massFluxFrac_PS(ccBurntSt.size(),0.);  

           forAll(ccBurntSt,faceI)
           {
             if (ccBurntSt[faceI] > 0.)  // unburnt val is -1.0 and burnt val i 1.0
             {
                if (UCBurntSt[faceI] > 0.)
                {
                 // Things to do when outer UC cell has burnt through
                 massFluxFrac_CC[faceI] = IUC_CCFrac; 
                 massFluxFrac_PS[faceI] = 1. - IUC_CCFrac; 
                 massFluxCCFracIUC_[localPyrolysisFaceI] = IUC_CCFrac; 
                }
                else
                {
                 // Things to do when outer unit cell is burning
                 if ( m_UC[faceI] > thresholdUCMass) 
                 {
                   massFluxFrac_CC[faceI] = 1.;
                   massFluxFrac_PS[faceI] = 0.;
                   massFluxCCFracOUC_[localPyrolysisFaceI] = 1.;
                 }
                 else
                 {
                   massFluxFrac_CC[faceI] = UC_CCFrac; 
                   massFluxFrac_PS[faceI] = 1. - UC_CCFrac; 
                   massFluxCCFracOUC_[localPyrolysisFaceI] = UC_CCFrac; 
                 }
                }
             }
             else
             {
                // nothing needs to be done, since massFluxFrac_CC by default is set to 1.
                //massFluxFrac_CC[faceI] = 1.;
                //massFluxFrac_PS[faceI] = 0.;
             }
             localPyrolysisFaceI++;
           }

           const label primaryPatchI = this->primaryPatchIDs()[i];

           fvPatchScalarField& fuelCC_patch = const_cast<fvPatchScalarField&>(fuelCC_bf[primaryPatchI]);
           fvPatchScalarField& fuelPS_patch = const_cast<fvPatchScalarField&>(fuelPS_bf[primaryPatchI]);

          if (isA<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelCC_patch) && 
             isA<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelPS_patch)  ) 
          {  
             Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField& 
                 nonUniFuelCC = refCast<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelCC_patch);

             nonUniFuelCC.setMassFluxFraction(massFluxFrac_CC); 

             Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField& 
                 nonUniFuelPS = refCast<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelPS_patch);

             nonUniFuelPS.setMassFluxFraction(massFluxFrac_PS); 
          }
          else
          {
        FatalErrorIn
        (
            "pyroCUPOneDimV1::updateFuelFluxMassFrac"
            "("
            ")"
        )
            << " Patch type not correct for patch \n"
            << gasMesh.boundaryMesh()[primaryPatchI].name() << nl
            << "for species " << speciesCC_ << " and " << speciesPS_ << "." << nl
            << " Please use patch type "
            << " nonUniFlowRateAdvectiveDiffusive "
            << "for both the species." << exit(FatalError);             
          }

    }
}


tmp<volScalarField> pyroCUPOneDimV1::kappaRad() const
{
    tmp<volScalarField> emiField(radiation_->absorptionEmission().e());
    return emiField;
}


void pyroCUPOneDimV1::setQrad(const scalarField& qrad, label patchI)
{
    qradBnd_.boundaryFieldRef()[patchI] = qrad;
}


void pyroCUPOneDimV1::setQconv(const scalarField& qconv, label patchI)
{
    qconvBnd_.boundaryFieldRef()[patchI] = qconv;
}


//void pyroCUPOneDimV1::solveEnergy()
//{
//    reactingOneDim21CharOxi::solveEnergy();
//}


void pyroCUPOneDimV1::preEvolveRegion()
{
    reactingOneDim21CharOxi::preEvolveRegion();

    const label vIndex = solidThermo_.composition().species()["v"];
    const label chIndex = solidThermo_.composition().species()["char"];
    const label dSpecIndex = solidThermo_.composition().species()["dSpec"];

    scalar rhoV = solidThermo_.composition().rho(vIndex,  10000, 298);
    scalar rhoC = solidThermo_.composition().rho(chIndex, 10000, 298);

     // Mass threshold to decide on transition from regime 1 to regime 2
     scalar thresholdUCMass = initMassUC_.value()-massFracUC_.value()*initMassUCCC_.value();

    label localPyrolysisFaceI = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];

        scalarField& ccBurntSt = swccGone_.boundaryFieldRef()[patchI];
        scalarField& UCBurntSt = outerUCGone_.boundaryFieldRef()[patchI];

        const scalarField& mUC = massUC_.boundaryField()[patchI];
        //scalarField& alpha_UC = alphaUC_.boundaryField()[patchI];           

        //const scalarField& cellV = regionMesh().V();

        const scalarField& mOLInit_patch   = mOLInit_.boundaryField()[patchI];
        const scalarField& mOLLost_patch   = mOLLost_.boundaryField()[patchI];

        scalarField& IOL_patch       = IOL_.boundaryFieldRef()[patchI];
        scalarField& IOUCHU_patch    = IOUCHU_.boundaryFieldRef()[patchI];
        scalarField& IOUCR1_patch    = IOUCR1_.boundaryFieldRef()[patchI];
        scalarField& IOUCR2_patch    = IOUCR2_.boundaryFieldRef()[patchI];
        scalarField& IIUC_patch      = IIUC_.boundaryFieldRef()[patchI];

        const scalarField& T_UC      = tempUC_.boundaryField()[patchI];           

        numPatchFaceCt_[i] = 0;
        facesOL_[i] = -1;
        facesOUC_[i] = -1;
        facesOUC2_[i] = -1;
        facesIUC_[i] = -1;

        scalar gasToVFac = rhoV/(rhoV-rhoC); //vapor pyrolsate to virgin mass conversion (i.e., 1 gm of vapor pyrolsate comes from this amount of virgin)..
        IOL_patch    = pos((1.-OLCrit_.value())*mOLInit_patch - gasToVFac*mOLLost_patch);

        ccBurntSt    = 1. - 2.*IOL_patch;

        IIUC_patch   = pos(0.01*initMassUC_.value() - mUC);

        UCBurntSt    = 2.*IIUC_patch - 1.;

        IOUCHU_patch = pos(ccBurntSt)*pos(igniTempUC_.value() - T_UC - SMALL); // *pos(-UCBurntSt);  // added pos(-UCBurntSt) later to get the test case to work

        scalarField thresCheck = pos(mUC - thresholdUCMass);
         
        IOUCR1_patch = pos(ccBurntSt)*pos(0.5 - IOUCHU_patch)*thresCheck; // pos(mUC - thresholdUCMass); 

        IOUCR2_patch = pos(-UCBurntSt)*pos(0.5 - thresCheck)*pos(0.5 - IOUCHU_patch); // pos(thresholdUCMass- mUC); 

        // If OUC-Regime 2 is active, then (based on remaining UC mass) compute energy split multiplier for OUC and IUC

        scalar OUC_Crit = 0.01*initMassUC_.value();
        scalar energySplit_thres = OUCEnSplit_.value()*thresholdUCMass;
        if (energySplit_thres < OUC_Crit)
        {
         energySplit_thres = OUC_Crit + 1.e-16; 
        }

        thresCheck = IOUCR2_patch; // used here as a dummy variable.. don't confuse with name..

        IOUCR2_patch = pos(thresCheck - 0.5)*max(pos(mUC - energySplit_thres),0.5);
        IIUC_patch = IIUC_patch + pos(thresCheck - 0.5)*(1. - IOUCR2_patch); 
       
        // -----------

       forAll(mUC,faceI)
       { 
         const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI];
         scalar IOLgone = 1-IOL_patch[faceI];
         forAll(cells,k)
         {
             Ys_[dSpecIndex][cells[k]] = IOLgone*1. + IOL_patch[faceI]*Ys_[dSpecIndex][cells[k]];
             Ys_[vIndex][cells[k]] = IOLgone*0. + IOL_patch[faceI]*Ys_[vIndex][cells[k]];
             Ys_[chIndex][cells[k]] = IOLgone*0. + IOL_patch[faceI]*Ys_[chIndex][cells[k]];
         }
         localPyrolysisFaceI = localPyrolysisFaceI + 1; 
       }


    }


}

// #####################################################
void pyroCUPOneDimV1::solveSpeciesMass()
{
    if (debug)
    {
        Info<< "pyroCUPOneDimV1::solveSpeciesMass()" << endl;
    }

    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho_, Yi)
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

// guptaa .. Adding a line below to limit the maximum value of Yi to 1
        Yi.min(1.0);
        Yt += Yi; 
    }

// guptaa .. Adding a line below to limit the maximum value of Yt to 1
    Yt.min(1.0);

    Ys_[Ys_.size() - 1] = 1.0 - Yt;

}

// #####################################################
void pyroCUPOneDimV1::evolveRegion()
{

    reactingOneDim21CharOxi::evolveRegion();

    //const label vIndex = solidThermo_.composition().species()["v"];
    //const label chIndex = solidThermo_.composition().species()["char"];
    //const label dSpecIndex = solidThermo_.composition().species()["dSpec"];

    //const scalarField& cellV = regionMesh().V();
    totalMassOL_  = 0.;

    // Mass threshold to decide on transition from regime 1 to regime 2
    //scalar thresholdUCMass = initMassUC_.value()-massFracUC_.value()*initMassUCCC_.value();

    scalarField& T = solidThermo_.T();

    volScalarField& T1 = solidThermo_.T();
    const volScalarField::Boundary& Tbf = T1.boundaryField();
  
    massFluxOL_  = 0.;
    massFluxOUC_ = 0.;
    massFluxIUC_ = 0.;
    massFluxOCC_ = 0.;
    massFluxOPS_ = 0.;
    massFluxICC_ = 0.;
    massFluxIPS_ = 0.;

    // Emissivity and absorptivity settings
    emmBnd_ = radiation_->absorptionEmission().e()();
    absBnd_ = emmBnd_;

    label localPyrolysisFaceI = 0;

    // Settings for updating fuel mass fraction
    const fvMesh& gasMesh = this->primaryMesh();
    const volScalarField& fuelCC = gasMesh.lookupObject<volScalarField>(speciesCC_); 
    const volScalarField::Boundary& fuelCC_bf = fuelCC.boundaryField();

    const volScalarField& fuelPS = gasMesh.lookupObject<volScalarField>(speciesPS_); 
    const volScalarField::Boundary& fuelPS_bf = fuelPS.boundaryField();

    massFluxCCFracIUC_ = 1.;
    massFluxCCFracOUC_ = 1.;

    scalar UC_CCFrac = (1.-massFracUC_.value())*initMassUCCC_.value()*hocPyrCC_.value() /
                   ( (1.-massFracUC_.value())*initMassUCCC_.value()*hocPyrCC_.value() + initMassUCPS_.value()*hocPyrPS_.value() );

    scalar IUC_CCFrac = initMassIUCCC_.value()*hocPyrCC_.value() /
                   ( initMassIUCCC_.value()*hocPyrCC_.value() + initMassIUCPS_.value()*hocPyrPS_.value() );

 
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i]; 

        const fvPatchScalarField& Tp = Tbf[patchI];

        // Getting the heat loss at the back boundary of the OL; to be used as heat source for OUC where OL is present
        const fvPatchScalarField& TbackBnd = Tbf[backBndPatch_[i]];
        const Foam::constHTemperatureFvPatchScalarField& TBCBackBnd = refCast<const Foam::constHTemperatureFvPatchScalarField>(TbackBnd);
        scalarField backBndHeatLoss = -TBCBackBnd.snGrad()*solidThermo_.kappa(backBndPatch_[i]);  // negavtive sign is used here because snGrad is negative..


        // Begin .. Getting scalar field for O2 in the adjacent gas-phase cell..
        // This is used to scale Q_flame contribution..

/*
        const fvPatch& patch = regionMesh().boundary()[patchI];

        // Get the coupling information from the mappedPatchBase
        const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        (
            patch.patch()
        );
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[mpp.samplePolyPatch().index()];

        fvPatchScalarField O2 = nbrPatch.lookupPatchField<volScalarField, scalar>("O2");
        scalarField alphaSgs  = nbrPatch.lookupPatchField<volScalarField, scalar>("alphaSgs");
        scalarField alpha     = nbrPatch.lookupPatchField<volScalarField, scalar>("thermo:alpha");    //for now use molecular alpha
        scalarField alphaDelta = (alpha+alphaSgs) * nbrPatch.deltaCoeffs();

        scalarField O2Int = O2.patchInternalField();
        mpp.distribute(O2Int);
        mpp.distribute(alphaDelta);

*/
        // End .. Getting scalar field for O2 in the adjacent gas-phase cell..


        const scalarField& qradPatch  = qradBnd_.boundaryField()[patchI];
        const scalarField& qconvPatch = qconvBnd_.boundaryField()[patchI];

        const scalarField& ccBurntSt = swccGone_.boundaryField()[patchI];
        //const scalarField& UCBurntSt = outerUCGone_.boundaryField()[patchI];

        scalarField& T_UC     = tempUC_.boundaryFieldRef()[patchI];           
        scalarField& m_UC     = massUC_.boundaryFieldRef()[patchI];           
        const scalarField& Cp_UC    = CpUC_.boundaryField()[patchI];           
        //scalarField& alpha_UC = alphaUC_.boundaryField()[patchI];           
        //scalarField& emm_UC   = emmUC_.boundaryField()[patchI];           

        scalarField& mdot_UC  = mdotUC_.boundaryFieldRef()[patchI];           
        scalarField& mdot_IUC = mdotIUC_.boundaryFieldRef()[patchI];           
        mdot_UC  = 0.;
        mdot_IUC = 0.;

        const scalarField& area_bFacet   = Tp.patch().magSf();

        scalarField& phiGasp = phiGas_.boundaryFieldRef()[patchI];

        scalarField& emmBndV   = emmBnd_.boundaryFieldRef()[patchI];           
        scalarField& absBndV   = absBnd_.boundaryFieldRef()[patchI];           

        scalarField massFluxFrac_CC(ccBurntSt.size(),1.); 
        scalarField massFluxFrac_PS(ccBurntSt.size(),0.);  

        scalarField qnet_wo_emm(ccBurntSt.size(),0.);  
        scalarField qemm(ccBurntSt.size(),0.);  
        scalarField qnet(ccBurntSt.size(),0.);  
        scalarField Tbnd(ccBurntSt.size(),298.);  
        scalarField deltaTemp(ccBurntSt.size(),0.);  

        //const scalarField& mOLInit_patch   = mOLInit_.boundaryField()[patchI];
        scalarField& mOLLost_patch         = mOLLost_.boundaryFieldRef()[patchI];

        const scalarField& IOL_patch       = IOL_.boundaryField()[patchI];
        const scalarField& IOUCHU_patch    = IOUCHU_.boundaryField()[patchI];
        const scalarField& IOUCR1_patch    = IOUCR1_.boundaryField()[patchI];
        const scalarField& IOUCR2_patch    = IOUCR2_.boundaryField()[patchI];
        const scalarField& IIUC_patch      = IIUC_.boundaryField()[patchI];

        scalarField& energyRelUC             = enerRelUC_.boundaryFieldRef()[patchI];

        //scalarField& pyrHOC_patch          = pyrolHOC_.boundaryFieldRef()[patchI];


        scalarField tempEmm = emmBndV;
        emmBndV = emmBndV*IOL_patch + emmIOUCHU_.value()*IOUCHU_patch + emmIOUCR1_.value()*IOUCR1_patch + emmIOUCR2_.value()*IOUCR2_patch + emmIUC_.value()*IIUC_patch; 

        // Coupling OL with OUC... Heat loss from OL back surface is given as heat source to the OUC
        T_UC = T_UC + IOL_patch*backBndHeatLoss*time_.deltaT().value()/(Cp_UC*(m_UC + 1.e-16));
        T_UC = max(298.,min(igniTempUC_.value()+5.,T_UC));     

        absBndV = absBndV*IOL_patch + absIOUCHU_.value()*IOUCHU_patch + absIOUCR1_.value()*IOUCR1_patch + absIOUCR2_.value()*IOUCR2_patch + absIUC_.value()*IIUC_patch; 
 
       scalarField QLoss(ccBurntSt.size(),0.); 
       scalarField QConvLoss(ccBurntSt.size(),0.); 
       scalarField effEmm(ccBurntSt.size(),0.); 

       scalarField qnet_to_IUC(ccBurntSt.size(),0.); 
       scalarField qnet_to_OUC(ccBurntSt.size(),0.); 
       scalarField qnet_to_OUCR1(ccBurntSt.size(),0.); 

       if (EnLossFracSpec_)
       {
         QLoss     = qradPatch*UCEnLossFr_.value();
       }
       else
       {
         QLoss     = UCEnLossFixed_.value();
       }


      if (ConvLossFracSpec_)
      {
         QConvLoss = qradPatch*(UCConvLossFr_.value()*IOUCR1_patch + UCConvLossFrR2_.value()*(IOUCR2_patch + IIUC_patch));
      }
      else
      {
         QConvLoss = UCConvLossFixed_.value();
      }
     
       qnet_wo_emm = (qradPatch*absBndV + 0.5*QLoss)*area_bFacet;

       scalarField IUCExist = pos(IIUC_patch - 0.25);
       scalarField OUCExist = pos(IOUCR2_patch - 0.25);

       qnet_to_OUCR1 = (qnet_wo_emm + (QFlameUC_.value()+QFlameExtraOUCR1_.value())*area_bFacet)*IOUCR1_patch;

       qnet_to_OUC = (qnet_wo_emm + (QFlameUC2_.value()+QFlameExtra_.value())*area_bFacet)*(1. - IUCExist) +  0.67*(qradPatch*absBndV + 0.2*QLoss + multFacQFl_.value()*QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet*IUCExist;
       //Line below is for CUP and the one after is for Class 3
       //qnet_to_OUC = (qnet_wo_emm + (QFlameUC2_.value()+QFlameExtra_.value())*area_bFacet)*(1. - IUCExist) +  0.67*(qradPatch*absBndV + 0.2*QLoss + 1.5*QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet*IUCExist;
       //qnet_to_OUC = (qnet_wo_emm + (QFlameUC2_.value()+QFlameExtra_.value())*area_bFacet)*(1. - IUCExist) +  0.67*(qradPatch*absBndV + 0.2*QLoss + QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet*IUCExist;
       qnet_to_OUC = qnet_to_OUC*OUCExist;

       qnet_to_IUC = (qnet_wo_emm + (1.5*QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet)*(1.-OUCExist) +  0.33*(qradPatch*absBndV + 0.4*QLoss + multFacQFl_.value()*QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet*OUCExist;
       //Line below is for CUP and the one after is for Class 3
       //qnet_to_IUC = (qnet_wo_emm + (1.5*QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet)*(1.-OUCExist) +  0.33*(qradPatch*absBndV + 0.4*QLoss + 1.5*QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet*OUCExist;
       //qnet_to_IUC = (qnet_wo_emm + (QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet)*(1.-OUCExist) +  0.33*(qradPatch*absBndV + 0.4*QLoss + QFlameUC2_.value() + QFlameExtra_.value())*area_bFacet*OUCExist;

       qnet_to_IUC = qnet_to_IUC*IUCExist;
 
       // Portion of incident radiation that is not consumed in the unit-cell
       scalarField qinc_rem(ccBurntSt.size(),0.);    
       scalar thresUCEmm = 25;  // threshold UC emission (i.e., minimum emission expected from UC) 
       qinc_rem = qradPatch - qradPatch*absBndV - QConvLoss - 0.5*QLoss -0.1*QLoss*IUCExist*OUCExist;

       // The above estimates are used to compute the energy being release inside the UC due to fuel burning
       scalarField QFlExtra = IOUCR1_patch*QFlameExtraOUCR1_.value() + QFlameExtra_.value()*(IOUCR2_patch + IIUC_patch); 
       scalarField QFlame   = IOUCR1_patch*QFlameUC_.value() + QFlameUC2_.value()*(IOUCR2_patch + IIUC_patch) + (multFacQFl_.value()-1.)*QFlameUC2_.value()*IUCExist; 
            // line below is for CUP and the one after is for Class 3
       //scalarField QFlame   = IOUCR1_patch*QFlameUC_.value() + QFlameUC2_.value()*(IOUCR2_patch + IIUC_patch) + 0.5*QFlameUC2_.value()*IUCExist; 
       //scalarField QFlame   = IOUCR1_patch*QFlameUC_.value() + QFlameUC2_.value()*(IOUCR2_patch + IIUC_patch); 
       energyRelUC = max(thresUCEmm + QFlame  - qinc_rem, QFlame ); 
     
 
       Tbnd  =  Tbnd*IOL_patch + T_UC*(IOUCHU_patch + IOUCR1_patch + IOUCR2_patch) +  tempIUC_.value()*IIUC_patch; 
       qemm  = emmBndV*Foam::constant::physicoChemical::sigma.value()*pow(Tbnd,4)*area_bFacet;       
       qnet  = (qradPatch*absBndV + qconvPatch)*area_bFacet - qemm;


       deltaTemp = qnet*time_.deltaT().value()/(Cp_UC*area_bFacet*(m_UC + 1.e-16));
       T_UC  = T_UC + deltaTemp*IOUCHU_patch;
 
       T_UC = max(298.,min(igniTempUC_.value()+5.,T_UC));     
       Tbnd  =  Tbnd*IOL_patch + T_UC*(IOUCHU_patch + IOUCR1_patch + IOUCR2_patch) +  tempIUC_.value()*IIUC_patch; 

       // Computing effective emissivity from emissive flux for UC and IUC.. 
       scalarField qemmUC(ccBurntSt.size(),0.);    
       qemmUC = qinc_rem + energyRelUC - QFlame;  
       effEmm = max(qemmUC,0.) / (Foam::constant::physicoChemical::sigma.value()*pow(Tbnd,4));

       emmBndV = emmBndV*IOL_patch + emmIOUCHU_.value()*IOUCHU_patch + effEmm*(IOUCR1_patch + IOUCR2_patch + IIUC_patch); 
       absBndV = absBndV*IOL_patch + absIOUCHU_.value()*IOUCHU_patch + 1.*(IOUCR1_patch + IOUCR2_patch + IIUC_patch); 

       
        scalar factorQ  =  1.; // min(alphaDelta*O2MassFrac/0.0024 , 1.); 
        scalar factorQ2 = 1.; 

       mdot_UC  = qnet_to_OUCR1/heatGassUC_.value() + qnet_to_OUC/heatGassUC2_.value();

       mdot_UC  = max(min(mdot_UC,m_UC*area_bFacet/time_.deltaT().value()),0.);
       m_UC     = max(m_UC - mdot_UC*time_.deltaT().value()/area_bFacet,0.);

       mdot_IUC = max(qnet_to_IUC/heatGassIUC_.value(),0.);       

       mOLLost_patch = mOLLost_patch + phiGasp*IOL_patch*time_.deltaT().value(); 

       forAll(m_UC,faceI)
       { 
         localPyrolysisFaceI = startFaceID_[i] + faceI;
         const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI];

         scalar TbndVal = Tbnd[faceI];
         scalar IOLgone = 1-IOL_patch[faceI];
         forAll(cells,k)
         {
             T[cells[k]] = IOLgone*TbndVal + T[cells[k]]*IOL_patch[faceI];
         }
         localPyrolysisFaceI = localPyrolysisFaceI + 1; 
       }


       if ( sum(mdot_IUC) > 1.e-16)
       {  
        forAll(massFluxIUC_,palletID)
        {
          scalarField& face2Pallet_patch   = face2Pallet_[palletID].boundaryFieldRef()[patchI];
          massFluxIUC_[palletID] = massFluxIUC_[palletID] + sum(mdot_IUC*face2Pallet_patch); 
        }
       }
       else
       {
         massFluxIUC_ = 0.;
       }

      mdot_UC   = max(mdot_UC - (energyRelUC)*area_bFacet*factorQ*IOUCR1_patch/hocPyrUC_.value() - (energyRelUC)*area_bFacet*factorQ2*IOUCR2_patch/hocPyrUC2_.value(),0.);

       phiGasp  = phiGasp*IOL_patch + mdot_UC + mdot_IUC;

    }


       scalarField mLRScaleRat(nPallets_,1.);
       scalarField totalM(nPallets_,0.);
       scalarField totalMLossRate(nPallets_,0.);
       scalarField totalPArea(nPallets_,0.);
       // getting global sum of the total pallet mass and the total mass loss rate
       totalM  = totalMassIUC_;
       combineReduce(totalM,sumOp<scalarField>());

       totalMLossRate  = massFluxIUC_;
       combineReduce(totalMLossRate,sumOp<scalarField>());

    if ( sum(totalMLossRate) > 1.e-16)
    {
       totalPArea = palletArea_;
       combineReduce(totalPArea,sumOp<scalarField>());

       forAll(totalMassIUC_,palletID)
       {
          if (totalMLossRate[palletID]*time_.deltaT().value() >totalM[palletID])
          {
              mLRScaleRat[palletID]  = totalM[palletID]/max(totalMLossRate[palletID]*time_.deltaT().value(),SMALL);
              massFluxIUC_[palletID] = massFluxIUC_[palletID] * mLRScaleRat[palletID];
          }
       }

    forAll(totalMassIUC_,palletID)
    {
       totalMassIUC_[palletID] = totalM[palletID] - mLRScaleRat[palletID]*totalMLossRate[palletID]*time_.deltaT().value();
    }

    localPyrolysisFaceI = 0;
    totalMassOUC_ = 0.;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI           = intCoupledPatchIDs_[i];
        //const scalarField& UCBurntSt = outerUCGone_.boundaryField()[patchI];
        scalarField& m_IUC           = massInnerUC_.boundaryFieldRef()[patchI];
        //const scalarField& mUC       = massUC_.boundaryField()[patchI];
        scalarField& phiGasp         = phiGas_.boundaryFieldRef()[patchI];
        scalarField& mdot_IUC        = mdotIUC_.boundaryFieldRef()[patchI];
        const scalarField& mdot_UC   = mdotUC_.boundaryField()[patchI];
        const scalarField& fArea  = regionMesh().boundary()[patchI].magSf();
        
   
        scalarField temp_field(m_IUC.size(),0.);
        m_IUC = 0.; 
           // Since inner unit cell mass tracking is based on global depletion, the mass
           // is updated even though the inner unit cell mass is not consumed at this patch face.
        forAll(massFluxIUC_,palletID)
        {
           scalarField& face2Pallet_patch   = face2Pallet_[palletID].boundaryFieldRef()[patchI];
           m_IUC = m_IUC + (totalMassIUC_[palletID]/totalPArea[palletID])*face2Pallet_patch;
           temp_field = temp_field + (mLRScaleRat[palletID])*face2Pallet_patch;
        }
        mdot_IUC = mdot_IUC*temp_field;
  

        const scalarField& IOL_patch       = IOL_.boundaryField()[patchI];
        //const scalarField& IOUCR1_patch    = IOUCR1_.boundaryField()[patchI];
        //const scalarField& IOUCR2_patch    = IOUCR2_.boundaryField()[patchI];
        const scalarField& IIUC_patch      = IIUC_.boundaryField()[patchI];

        const scalarField& energyRelUC     = enerRelUC_.boundaryField()[patchI];

        mdot_IUC   = max(mdot_IUC - (energyRelUC)*fArea*IIUC_patch/hocPyrIUC_.value(),0.);

        phiGasp  = phiGasp*IOL_patch + mdot_UC  + mdot_IUC;
             

    }

  }


    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI                 = intCoupledPatchIDs_[i];
        const scalarField& mdot_IUC        = mdotIUC_.boundaryField()[patchI];
        const scalarField& mdot_UC         = mdotUC_.boundaryField()[patchI];
        const scalarField& IOL_patch       = IOL_.boundaryField()[patchI];
        const scalarField& IOUCHU_patch    = IOUCHU_.boundaryField()[patchI];
        const scalarField& IOUCR1_patch    = IOUCR1_.boundaryField()[patchI];
        const scalarField& IOUCR2_patch    = IOUCR2_.boundaryField()[patchI];
        const scalarField& IIUC_patch      = IIUC_.boundaryField()[patchI];
        scalarField& pyrHOC_patch          = pyrolHOC_.boundaryFieldRef()[patchI];

        scalarField massFluxFrac_CC(IOL_patch.size(),1.); 
        scalarField massFluxFrac_PS(IOL_patch.size(),0.);  

       pyrHOC_patch = hocPyrUC_.value()*IOUCR1_patch + (IOUCR2_patch + IIUC_patch)*(mdot_UC*hocPyrUC2_.value() + mdot_IUC*hocPyrIUC_.value())/(mdot_UC + mdot_IUC + 1.e-16);  

       massFluxFrac_CC = 1.*IOL_patch + 1.*IOUCHU_patch + 1.*IOUCR1_patch + (IOUCR2_patch + IIUC_patch)*(mdot_UC*hocPyrUC2_.value()*UC_CCFrac + mdot_IUC*hocPyrIUC_.value()*IUC_CCFrac)/(mdot_UC*hocPyrUC2_.value() + mdot_IUC*hocPyrIUC_.value() + 1.e-16);

       massFluxFrac_PS = 1. - massFluxFrac_CC;

           // Calling function for the gas-phase species BC to set fuel mass fractions for incoming pyrolysate..
           const label primaryPatchI = this->primaryPatchIDs()[i];

           fvPatchScalarField& fuelCC_patch = const_cast<fvPatchScalarField&>(fuelCC_bf[primaryPatchI]);
           fvPatchScalarField& fuelPS_patch = const_cast<fvPatchScalarField&>(fuelPS_bf[primaryPatchI]);

          if (isA<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelCC_patch) &&
             isA<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelPS_patch)  )
          {  
             Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField& 
                 nonUniFuelCC = refCast<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelCC_patch);

             nonUniFuelCC.setMassFluxFraction(massFluxFrac_CC); 

             Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField&  
                nonUniFuelPS = refCast<Foam::nonUniFlowRateAdvectiveDiffusiveFvPatchScalarField>(fuelPS_patch);

             nonUniFuelPS.setMassFluxFraction(massFluxFrac_PS); 
          }
          else
          {
        FatalErrorIn
        (
            "pyroCUPOneDimV1::updateFuelFluxMassFrac"
            "("
            ")"
        )
            << " Patch type not correct for patch \n"
            << gasMesh.boundaryMesh()[primaryPatchI].name() << nl
            << "for species " << speciesCC_ << " and " << speciesPS_ << "." << nl
            << " Please use patch type "            << " nonUniFlowRateAdvectiveDiffusive "
            << "for both the species." << exit(FatalError);             
          }

           // Calling function for the back boundary of OL to update the Tref  
           const fvPatchScalarField& TbackBnd = Tbf[backBndPatch_[i]];
           const Foam::constHTemperatureFvPatchScalarField& TBCBackBnd = refCast<const Foam::constHTemperatureFvPatchScalarField>(TbackBnd);
           const scalarField& TUC = tempUC_.boundaryField()[patchI];
           constHTemperatureFvPatchScalarField& TBCBack = const_cast<constHTemperatureFvPatchScalarField&>(TBCBackBnd);
           TBCBack.setTInf(TUC);

    }

    // Reduce operations to get the mass loss rates for pallets

    reduce(massFluxOL_,sumOp<scalarField>());
    reduce(massFluxOUC_,sumOp<scalarField>());
    reduce(massFluxIUC_,sumOp<scalarField>());
    reduce(massFluxOCC_,sumOp<scalarField>());
    reduce(massFluxOPS_,sumOp<scalarField>());
    reduce(massFluxICC_,sumOp<scalarField>());
    reduce(massFluxIPS_,sumOp<scalarField>());

    // Reduce operations to get the total remaining mass of pallets

    diagTMOL_  = totalMassOL_;
    diagTMOUC_ = totalMassOUC_;
    diagTMIUC_ = totalMassIUC_;

    reduce(diagTMOL_,sumOp<scalarField>());
    reduce(diagTMOUC_,sumOp<scalarField>());
    reduce(diagTMIUC_,sumOp<scalarField>());

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pyroCUPOneDimV1::pyroCUPOneDimV1
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    reactingOneDim21CharOxi(modelType, mesh, regionType),

    tempUC_
    (
        IOobject
        (
            "tempUC",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    massUC_
    (
        IOobject
        (
            "massUC",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    massUCCC_
    (
        IOobject
        (
            "massUCCC",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        massUC_
    ),

    massUCPS_
    (
        IOobject
        (
            "massUCPS",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        massUC_
    ),

    CpUC_
    (
        IOobject
        (
            "CpUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimMass/dimTemperature, 0.0)
    ),

    alphaUC_
    (
        IOobject
        (
            "alphaUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    emmUC_
    (
        IOobject
        (
            "emmUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mdotUC_
    (
        IOobject
        (
            "mdotUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    mdotIUC_
    (
        IOobject
        (
            "mdotIUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    swccGone_
    (
        IOobject
        (
            "swccGone",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, -1.0)
    ),

    outerUCGone_
    (
        IOobject
        (
            "outerUCGone",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, -1.0)
    ),

    massInnerUC_
    (
        IOobject
        (
            "massInnerUC",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    heatGassUC_(dimensionedScalar("heatGassUC", dimEnergy/dimMass, 0.0)),
    heatGassUC2_(dimensionedScalar("heatGassUC2", dimEnergy/dimMass, 0.0)),
    QFlameUC_(dimensionedScalar("QFlameUC", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    QFlameUC2_(dimensionedScalar("QFlameUC2", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    QEmmUC_(dimensionedScalar("QEmmUC", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    OLCrit_(dimensionedScalar("OLCrit", dimless, 0.0)),
    OUCEnSplit_(dimensionedScalar("OUCEnSplit", dimless, 0.02)),
    multFacQFl_(dimensionedScalar("multFacQFl", dimless, 1.)),
    massFracUC_(dimensionedScalar("CCMassFracUC", dimless, 0.6)),
    TwoHeatGass_(false),
    Tacti_(dimensionedScalar("Tacti", dimTemperature, 0.0)),
    AGassReac_(dimensionedScalar("AGassReac", dimMass/dimTime, 0.0)),
    nGassReac_(dimensionedScalar("nGassReac", dimless, 0.0)),
    initMassUC_(dimensionedScalar("initMassUC", dimMass/dimLength/dimLength, 0.0)),
    initMassUCCC_(dimensionedScalar("initMassUCCC", dimMass/dimLength/dimLength, 0.0)),
    initMassUCPS_(dimensionedScalar("initMassUCPS", dimMass/dimLength/dimLength, 0.0)),

    heatGassIUC_(dimensionedScalar("heatGassIUC", dimEnergy/dimMass, 0.0)),
    tempIUC_(dimensionedScalar("tempIUC", dimTemperature, 0.0)),
    emmIUC_(dimensionedScalar("emmIUC", dimless, 0.0)),
    absIUC_(dimensionedScalar("absIUC", dimless, 0.0)),
    initMassIUC_(dimensionedScalar("initMassIUC", dimMass/dimLength/dimLength, 0.0)),
    initMassIUCCC_(dimensionedScalar("initMassIUCCC", dimMass/dimLength/dimLength, 0.0)),
    initMassIUCPS_(dimensionedScalar("initMassIUCPS", dimMass/dimLength/dimLength, 0.0)),
    massIUC_(dimensionedScalar("massIUC", dimMass/dimLength/dimLength, 0.0)),

    igniTempUC_(dimensionedScalar("igniTempUC", dimTemperature, 0.0)),
    pModelTypeUC_(1),

    hocPyrUC_(dimensionedScalar("hocPyrUC", dimEnergy/dimMass, 0.0)),
    hocPyrUC2_(dimensionedScalar("hocPyrUC2", dimEnergy/dimMass, 0.0)),
    hocPyrIUC_(dimensionedScalar("hocPyrIUC", dimEnergy/dimMass, 0.0)),

    qradBnd_
    (
        IOobject
        (
            "qradBnd",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimLength/dimLength, 0.0)
    ),

    qconvBnd_
    (
        IOobject
        (
            "qconvBnd",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimLength/dimLength, 0.0)
    ),

   multiFuel_(false),
   speciesCC_("none"),
   speciesPS_("none"),

   nPallets_(8),

   nXPallets_(2),
   nYPallets_(2),
   nZPallets_(2),

   dXPallets_(dimensionedScalar("dXPallets", dimLength, 1.2192)),
   dYPallets_(dimensionedScalar("dYPallets", dimLength, 1.2192)),
   dZPallets_(dimensionedScalar("dZPallets", dimLength, 1.524)),

   lXPallet_(dimensionedScalar("lXPallet", dimLength, 1.0668)),
   lYPallet_(dimensionedScalar("lYPallet", dimLength, 1.0668)),
   lZPallet_(dimensionedScalar("lZPallet", dimLength, 1.0668)),

   botCorPalletOrigin_(dimensionedVector("botCorPalletOrigin", dimLength, Foam::vector(0,0,0))),

    IOL_
    (
        IOobject
        (
            "IOL",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IOUCHU_
    (
        IOobject
        (
            "IOUCHU",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IOUCR1_
    (
        IOobject
        (
            "IOUCR1",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IOUCR2_
    (
        IOobject
        (
            "IOUCR2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IIUC_
    (
        IOobject
        (
            "IIUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mOLInit_
    (
        IOobject
        (
            "mOLInit",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mOLLost_
    (
        IOobject
        (
            "mOLLost",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),


    emmIOUCHU_(dimensionedScalar("emmIOUCHU", dimless, 0.0)),
    absIOUCHU_(dimensionedScalar("absIOUCHU", dimless, 0.0)),
    emmIOUCR1_(dimensionedScalar("emmIOUCR1", dimless, 0.0)),
    absIOUCR1_(dimensionedScalar("absIOUCR1", dimless, 0.0)),
    emmIOUCR2_(dimensionedScalar("emmIOUCR1", dimless, 0.0)),
    absIOUCR2_(dimensionedScalar("absIOUCR2", dimless, 0.0)),
 
    UCEnLossFr_(dimensionedScalar("UCEnLossFr", dimless, 0.2)),
    EnLossFracSpec_(true),
    UCEnLossFixed_(dimensionedScalar("UCEnLossFixed", dimEnergy/dimTime/dimLength/dimLength, 20000.)),
    UCConvLossFr_(dimensionedScalar("UCConvLossFr", dimless, 0.1)),
    UCConvLossFrR2_(dimensionedScalar("UCConvLossFrR2", dimless, 0.1)),
    ConvLossFracSpec_(true),
    UCConvLossFixed_(dimensionedScalar("UCConvLossFixed", dimEnergy/dimTime/dimLength/dimLength, 10000.)),

    QFlameExtra_(dimensionedScalar("QFlameExtra", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    QFlameExtraOUCR1_(dimensionedScalar("QFlameExtraOUCR1", dimEnergy/dimTime/dimLength/dimLength, 0.0)),

    enerRelUC_
    (
        IOobject
        (
            "enerRelUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimLength/dimLength, 0.0)
    ),

    pyrolHOC_
    (
        IOobject
        (
            "pyrolHOC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
    )


{
   if (active_)
    {
        initParams();
    }
}


pyroCUPOneDimV1::pyroCUPOneDimV1
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    reactingOneDim21CharOxi(modelType, mesh, dict, regionType),

    tempUC_
    (
        IOobject
        (
            "tempUC",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    massUC_
    (
        IOobject
        (
            "massUC",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    massUCCC_
    (
        IOobject
        (
            "massUCCC",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        massUC_
    ),

    massUCPS_
    (
        IOobject
        (
            "massUCPS",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        massUC_
    ),

    CpUC_
    (
        IOobject
        (
            "CpUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimMass/dimTemperature, 0.0)
    ),

    alphaUC_
    (
        IOobject
        (
            "alphaUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    emmUC_
    (
        IOobject
        (
            "emmUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mdotUC_
    (
        IOobject
        (
            "mdotUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    mdotIUC_
    (
        IOobject
        (
            "mdotIUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    swccGone_
    (
        IOobject
        (
            "swccGone",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, -1.0)
    ),

    outerUCGone_
    (
        IOobject
        (
            "outerUCGone",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, -1.0)
    ),

    massInnerUC_
    (
        IOobject
        (
            "massInnerUC",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    heatGassUC_(dimensionedScalar("heatGassUC", dimEnergy/dimMass, 0.0)),
    heatGassUC2_(dimensionedScalar("heatGassUC2", dimEnergy/dimMass, 0.0)),
    QFlameUC_(dimensionedScalar("QFlameUC", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    QFlameUC2_(dimensionedScalar("QFlameUC2", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    QEmmUC_(dimensionedScalar("QEmmUC", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    OLCrit_(dimensionedScalar("OLCrit", dimless, 0.0)),
    OUCEnSplit_(dimensionedScalar("OUCEnSplit", dimless, 0.02)),
    multFacQFl_(dimensionedScalar("multFacQFl", dimless, 1.)),
    massFracUC_(dimensionedScalar("CCMassFracUC", dimless, 0.6)),
    TwoHeatGass_(false),
    Tacti_(dimensionedScalar("Tacti", dimTemperature, 0.0)),
    AGassReac_(dimensionedScalar("AGassReac", dimMass/dimTime, 0.0)),
    nGassReac_(dimensionedScalar("nGassReac", dimless, 0.0)),
    initMassUC_(dimensionedScalar("initMassUC", dimMass/dimLength/dimLength, 0.0)),
    initMassUCCC_(dimensionedScalar("initMassUCCC", dimMass/dimLength/dimLength, 0.0)),
    initMassUCPS_(dimensionedScalar("initMassUCPS", dimMass/dimLength/dimLength, 0.0)),

    heatGassIUC_(dimensionedScalar("heatGassIUC", dimEnergy/dimMass, 0.0)),
    tempIUC_(dimensionedScalar("tempIUC", dimTemperature, 0.0)),
    emmIUC_(dimensionedScalar("emmIUC", dimless, 0.0)),
    absIUC_(dimensionedScalar("absIUC", dimless, 0.0)),
    initMassIUC_(dimensionedScalar("initMassIUC", dimMass/dimLength/dimLength, 0.0)),
    initMassIUCCC_(dimensionedScalar("initMassIUCCC", dimMass/dimLength/dimLength, 0.0)),
    initMassIUCPS_(dimensionedScalar("initMassIUCPS", dimMass/dimLength/dimLength, 0.0)),
    massIUC_(dimensionedScalar("massIUC", dimMass/dimLength/dimLength, 0.0)),

    igniTempUC_(dimensionedScalar("igniTempUC", dimTemperature, 0.0)),
    pModelTypeUC_(1),
    hocPyrUC_(dimensionedScalar("hocPyrUC", dimEnergy/dimMass, 0.0)),
    hocPyrUC2_(dimensionedScalar("hocPyrUC2", dimEnergy/dimMass, 0.0)),
    hocPyrIUC_(dimensionedScalar("hocPyrIUC", dimEnergy/dimMass, 0.0)),

    qradBnd_
    (
        IOobject
        (
            "qradBnd",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimLength/dimLength, 0.0)
    ),

    qconvBnd_
    (
        IOobject
        (
            "qconvBnd",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimLength/dimLength, 0.0)
    ),

   multiFuel_(false),
   speciesCC_("none"),
   speciesPS_("none"),

   nPallets_(8),

   nXPallets_(2),
   nYPallets_(2),
   nZPallets_(2),

   dXPallets_(dimensionedScalar("dXPallets", dimLength, 1.2192)),
   dYPallets_(dimensionedScalar("dYPallets", dimLength, 1.2192)),
   dZPallets_(dimensionedScalar("dZPallets", dimLength, 1.524)),

   lXPallet_(dimensionedScalar("lXPallet", dimLength, 1.0668)),
   lYPallet_(dimensionedScalar("lYPallet", dimLength, 1.0668)),
   lZPallet_(dimensionedScalar("lZPallet", dimLength, 1.0668)),

   botCorPalletOrigin_(dimensionedVector("botCorPalletOrigin", dimLength, Foam::vector(0,0,0))),

    IOL_
    (
        IOobject
        (
            "IOL",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IOUCHU_
    (
        IOobject
        (
            "IOUCHU",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IOUCR1_
    (
        IOobject
        (
            "IOUCR1",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IOUCR2_
    (
        IOobject
        (
            "IOUCR2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    IIUC_
    (
        IOobject
        (
            "IIUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mOLInit_
    (
        IOobject
        (
            "mOLInit",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mOLLost_
    (
        IOobject
        (
            "mOLLost",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    emmIOUCHU_(dimensionedScalar("emmIOUCHU", dimless, 0.0)),
    absIOUCHU_(dimensionedScalar("absIOUCHU", dimless, 0.0)),
    emmIOUCR1_(dimensionedScalar("emmIOUCR1", dimless, 0.0)),
    absIOUCR1_(dimensionedScalar("absIOUCR1", dimless, 0.0)),
    emmIOUCR2_(dimensionedScalar("emmIOUCR1", dimless, 0.0)),
    absIOUCR2_(dimensionedScalar("absIOUCR2", dimless, 0.0)),

    UCEnLossFr_(dimensionedScalar("UCEnLossFr", dimless, 0.2)),
    EnLossFracSpec_(true),
    UCEnLossFixed_(dimensionedScalar("UCEnLossFixed", dimEnergy/dimTime/dimLength/dimLength, 20000.)),
    UCConvLossFr_(dimensionedScalar("UCConvLossFr", dimless, 0.1)),
    UCConvLossFrR2_(dimensionedScalar("UCConvLossFrR2", dimless, 0.1)),
    ConvLossFracSpec_(true),
    UCConvLossFixed_(dimensionedScalar("UCConvLossFixed", dimEnergy/dimTime/dimLength/dimLength, 10000.)),

    QFlameExtra_(dimensionedScalar("QFlameExtra", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    QFlameExtraOUCR1_(dimensionedScalar("QFlameExtraOUCR1", dimEnergy/dimTime/dimLength/dimLength, 0.0)),
    
    enerRelUC_
    (
        IOobject
        (
            "enerRelUC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimLength/dimLength, 0.0)
    ),

    pyrolHOC_
    (
        IOobject
        (
            "pyrolHOC",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
    )


{
   if (active_)
    {
        initParams();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pyroCUPOneDimV1::~pyroCUPOneDimV1()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void pyroCUPOneDimV1::info() const
{
    Info<< "\nCUP Pyrolysis in region: " << regionMesh().name() << endl;

    Info<< "Outer-liner MLR [kg/s]: " << time_.timeOutputValue();
    forAll(massFluxOL_,pID)
    {
        Info << " " << massFluxOL_[pID];
    }
    Info << nl;

    Info<< "Outer-UC MLR [kg/s]: " << time_.timeOutputValue();
    forAll(massFluxOUC_,pID)
    {
        Info << " " << massFluxOUC_[pID];
    }
    Info << nl;

    Info<< "Outer-UC-CC MLR [kg/s]: " << time_.timeOutputValue();
    forAll(massFluxOCC_,pID)
    {
        Info << " " << massFluxOCC_[pID];
    }
    Info << nl;

    Info<< "Outer-UC-PS MLR [kg/s]: " << time_.timeOutputValue();
    forAll(massFluxOPS_,pID)
    {
        Info << " " << massFluxOPS_[pID];
    }
    Info << nl;

    Info<< "Inner-UC MLR [kg/s]: " << time_.timeOutputValue();
    forAll(massFluxIUC_,pID)
    {
        Info << " " << massFluxIUC_[pID];
    }
    Info << nl;

    Info<< "Inner-UC-CC MLR [kg/s]: " << time_.timeOutputValue();
    forAll(massFluxICC_,pID)
    {
        Info << " " << massFluxICC_[pID];
    }
    Info << nl;

    Info<< "Inner-UC-PS MLR [kg/s]: " << time_.timeOutputValue();
    forAll(massFluxIPS_,pID)
    {
        Info << " " << massFluxIPS_[pID];
    }
    Info << nl;

    Info<< "Outer-Liner Mass [kg]: " << time_.timeOutputValue();
    forAll(diagTMOL_,pID)
    {
        Info << " " << diagTMOL_[pID];
    }
    Info << nl;

    Info<< "Outer UC Mass [kg]: " << time_.timeOutputValue();
    forAll(diagTMOUC_,pID)
    {
        Info << " " << diagTMOUC_[pID];
    }
    Info << nl;

    Info<< "Inner UC Mass [kg]: " << time_.timeOutputValue();
    forAll(diagTMIUC_,pID)
    {
        Info << " " << diagTMIUC_[pID];
    }
    Info << nl;

    /*
    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
     */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //

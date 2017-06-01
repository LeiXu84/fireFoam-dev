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

\*---------------------------------------------------------------------------*/

#include "physicoChemicalConstants.H"
#include "fvDOM.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ThermoCloud<CloudType>&
Foam::ThermoCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ThermoCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ThermoCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::SLGThermo& Foam::ThermoCloud<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ThermoCloud<CloudType>::T() const
{
    return T_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ThermoCloud<CloudType>::p() const
{
    return p_;
}


template<class CloudType>
inline const Foam::HeatTransferModel<Foam::ThermoCloud<CloudType>>&
Foam::ThermoCloud<CloudType>::heatTransfer() const
{
    return heatTransferModel_;
}


template<class CloudType>
inline const Foam::scalarIntegrationScheme&
Foam::ThermoCloud<CloudType>::TIntegrator() const
{
    return TIntegrator_;
}


template<class CloudType>
inline bool Foam::ThermoCloud<CloudType>::radiation() const
{
    return radiation_;
}

template<class CloudType>
inline bool Foam::ThermoCloud<CloudType>::coupledRadiation() const
{
    return coupledRadiation_;
}

template<class CloudType>
inline Foam::label Foam::ThermoCloud<CloudType>::nBands() const
{
    return nBands_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaP(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaP_[bandI];
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaP(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaP_[bandI];
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPSc(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPSc_[bandI];
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPSc(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPSc_[bandI];
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radT4(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radT4_[bandI];
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radT4(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radT4_[bandI];
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPT4(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPT4_[bandI];
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPT4(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPT4_[bandI];
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPScAsy(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPScAsy_[bandI];
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPScAsy(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPScAsy_[bandI];
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::CScat(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return CScat_[bandI];
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::CScat(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return CScat_[bandI];
}

template<class CloudType>
inline const Foam::List<Foam::scalar>& Foam::ThermoCloud<CloudType>::diaVal() const
{
  return diaVal_;
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsTrans()
{
    return hsTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsTrans() const
{
    return hsTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsCoeff()
{
    return hsCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsCoeff() const
{
    return hsCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermoCloud<CloudType>::Sh(volScalarField& hs) const
{
    if (debug)
    {
        Info<< "hsTrans min/max = " << min(hsTrans()).value() << ", "
            << max(hsTrans()).value() << nl
            << "hsCoeff min/max = " << min(hsCoeff()).value() << ", "
            << max(hsCoeff()).value() << endl;
    }

    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("h"))
        {
            const volScalarField Cp(thermo_.thermo().Cp());
            const volScalarField::Internal
                Vdt(this->mesh().V()*this->db().time().deltaT());

            return
                hsTrans()/Vdt
              - fvm::SuSp(hsCoeff()/(Cp*Vdt), hs)
              + hsCoeff()/(Cp*Vdt)*hs;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = -hsTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(hs, dimEnergy/dimTime));
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::Ep(const label bandI) const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    if (radiation_)
    {
        scalarField& Ep = tEp.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();
        const scalarField& sumAreaPT4 = radAreaPT4_[bandI].field();

        //Ep = sumAreaPT4*epsilon*physicoChemical::sigma.value()/V/dt;
        Ep = sumAreaPT4*physicoChemical::sigma.value()/dt;    // In this implementation, epsilon (i.e., emission efficiency) and division by V is included in radAreaP_ term
    }
    //Pout << "Sum of Ep: " << Foam::sum(Ep) << endl;
    return tEp;
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::ap(const label bandI) const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    if (radiation_)
    {
        scalarField& ap = tap.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& V = this->mesh().V();
        //const scalar epsilon = constProps_.epsilon0();
        const scalarField& sumAreaP = radAreaP_[bandI].field();

        //ap = sumAreaP*epsilon/V/dt;
        ap = sumAreaP/dt;     // In this implementation, epsilon (i.e., absorption efficiency) and division by V is included in radAreaP_ term
    }

    return tap;
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<CloudType>::sigmap(const label bandI) const
{
    tmp<volScalarField> tsigmap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:sigmap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    if (radiation_)
    {
        scalarField& sigmap = tsigmap.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();
        const scalar f = constProps_.f0();
        const scalarField& sumAreaP = radAreaPSc_[bandI].field();

        //sigmap *= sumAreaP*(1.0 - f)*(1.0 - epsilon)/V/dt;
        sigmap = sumAreaP/dt;     // In this implementation, scattering efficiency and division by V is included in radAreaP_ term
    }

    return tsigmap;
}

template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::phaseFunc(const label bandI, const label sour, const label dest) const
{
    tmp<volScalarField> tPhaseFunc
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:phaseFunc",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    if ((radiation_) && (sour != dest))
    {
        scalarField& phaseF = tPhaseFunc.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& cscat = CScat_[bandI].field();
        scalarField asyFacVal = radAreaPScAsy_[bandI].field()/(radAreaPSc_[bandI].field() + VSMALL);    
        scalarField A_asy = 2.*asyFacVal/(1. + asyFacVal*asyFacVal);
        scalarField B_asy = (1. - asyFacVal*asyFacVal)/(1. + asyFacVal*asyFacVal);

        const Foam::radiation::radiationModel& radModel_ = this->mesh().template lookupObject<radiation::radiationModel>("radiationProperties");
        const Foam::radiation::fvDOM& fvdomModel_=refCast<const Foam::radiation::fvDOM>(radModel_);
        const Foam::vector dRayS = fvdomModel_.IRay(sour).d();
        const Foam::vector dRayD = fvdomModel_.IRay(dest).d();
        scalar sDotd = dRayS & dRayD;

        //Info << "cscat: " << cscat << endl;
        phaseF = cscat * B_asy * (1. + A_asy*sDotd + A_asy*A_asy*sDotd*sDotd + Foam::pow(A_asy,3.)*Foam::pow(sDotd,3.) + Foam::pow(A_asy,4.)*Foam::pow(sDotd,4.)) / (4.*constant::mathematical::pi - fvdomModel_.IRay(sour).omega() ); 
     }

    return tPhaseFunc;
}

template<class CloudType>
inline void Foam::ThermoCloud<CloudType>::updateRadProperties(const label celli, const scalar ap, const label np, const scalar T4, const scalar dt)
{
        scalar pDia = pow((4.*ap/3.1415),0.5);
        forAll(radAreaP_, bandI)
        {       
		radAreaP_[bandI].field()[celli] += dt*np*ap*getQAbs(pDia,bandI);
		radT4_[bandI].field()[celli] += dt*np*T4;
		radAreaPT4_[bandI].field()[celli] += dt*np*ap*T4;
        }
}


template<class CloudType>
inline Foam::List<Foam::scalar> Foam::ThermoCloud<CloudType>::radCalc(const label celli, const scalar ap, const scalar np, const scalar T4, const scalar dt)
{
        const Foam::radiation::radiationModel& radModel_ = this->mesh().template lookupObject<radiation::radiationModel>("radiationProperties");
        const Foam::radiation::fvDOM& fvdomModel_=refCast<const Foam::radiation::fvDOM>(radModel_);        

        List<scalar> kgAndEmm(2,0.);
        scalar pDia = pow((4.*ap/3.1415),0.5);

        List<scalar> QAbsorp(nBands_,0.);
        List<scalar> QScat(nBands_,0.);
        List<scalar> asyFactor(nBands_,0.);

        //Info << " nBands: " << nBands_ << endl;
        //Info << "diaVal: " << diaVal() << endl;
        //Info << "radProp_ : " << radProp_ << endl;
        //Info << "numPropBands: " << numPropBands_ << endl;

        //Pout << "np celli pDia: " << np <<" " << celli << " " << pDia << endl;

        scalar npDen = np/this->mesh().V()[celli]; 

         //Pout << "npDen: " << npDen << endl;   
         //Pout << "QAbsorp: " << QAbsorp << endl;
         //Pout << "QScat: " << QScat << endl;
         //Pout << "asyFactor: " << asyFactor << endl;

        getQAbsQSct(pDia,npDen,QAbsorp,QScat,asyFactor); 
 
        //Pout << "here: " << " 1 " << endl;

        List<scalar> QEmm(nBands_,0.);
        List<scalar> emFrac(nBands_,0.);
        scalar pT = pow(T4,0.5); 
        getQEmmFrac(pDia,pT,QEmm,emFrac);   //NEED TO IMPLEMENT THIS FUNCTION....
    
        scalar multFac = dt*npDen*ap;
        forAll(radAreaP_, bandI)
        {      
                scalar QAbs = QAbsorp[bandI]; 
                scalar QSct = QScat[bandI]; 
                scalar asy  = asyFactor[bandI]; 
                scalar effEmm = QEmm[bandI]*emFrac[bandI];
		radAreaP_[bandI].field()[celli] += multFac*QAbs;
		radAreaPSc_[bandI].field()[celli] += multFac*QSct;
                //Info << " Abs contribution: " << dt*npDen*ap*QAbs  << endl;
		radT4_[bandI].field()[celli] += dt*npDen*T4;
		radAreaPT4_[bandI].field()[celli] += multFac*effEmm*T4;
		radAreaPScAsy_[bandI].field()[celli] += multFac*QSct*asy;

                kgAndEmm[0] += QAbs*fvdomModel_.GLambda(bandI)[celli];
                kgAndEmm[1] += effEmm;
 
        }
     
//        scalar QAbs = QAbsorp[0];
//        scalar absTemp = ap*kgAndEmm[0]*np*dt*1.e6;
//        scalar emmTemp = physicoChemical::sigma.value()*dt*np*ap*QAbs*T4*1.e6*fvdomModel_.GLambda(0)[celli]/(fvdomModel_.GLambda(0)[celli] + VSMALL);
//        Pout << "Abs by droplets here: " << absTemp << " " << emmTemp << endl;               

        //kgAndEmm[1] = getEmm(pDia);
        return kgAndEmm;     
}

template<class CloudType>
inline Foam::scalar Foam::ThermoCloud<CloudType>::getkG(const scalar pDia) const
{
scalar test1;

test1 = 0;

return test1;
} 

template<class CloudType>
inline void Foam::ThermoCloud<CloudType>::getQEmmFrac(const scalar pDia, const scalar pT, List<scalar>& QEmm, List<scalar>& emFrac) const
{
// NEED TO DEFINE THIS FUNCTION........
return;
}

template<class CloudType>
inline Foam::scalar Foam::ThermoCloud<CloudType>::getQAbs(const scalar pDia, const label bandI) const
{

scalar test1;

test1 = 0;

return test1;
}

template<class CloudType>
inline void Foam::ThermoCloud<CloudType>::getQAbsQSct(const scalar pDia, const scalar np, List<scalar>& QAbs, List<scalar>& QSct, List<scalar>& asyFac) const
{

    label begin = 0;
    label end = diaVal_.size()-1;
    label nDias = end + 1;
     
    List<scalar> qabs(numPropBands_,0.);
    List<scalar> qsct(numPropBands_,0.);
    List<scalar> asyF(numPropBands_,0.);

    if(nDias==1)
    {
        qabs = absEff_[begin];
        qsct = sctEff_[begin];
        asyF = asyFac_[begin];
    }
    else if(pDia <= diaVal_[begin])
    {
        qabs = (absEff_[begin+1]-absEff_[begin])/(diaVal_[begin+1]-diaVal_[begin])*(pDia - diaVal_[begin])+absEff_[begin];
        qsct = (sctEff_[begin+1]-sctEff_[begin])/(diaVal_[begin+1]-diaVal_[begin])*(pDia - diaVal_[begin])+sctEff_[begin];
        asyF = (asyFac_[begin+1]-asyFac_[begin])/(diaVal_[begin+1]-diaVal_[begin])*(pDia - diaVal_[begin])+asyFac_[begin];
    }
    else if(pDia >= diaVal_[end])
    {
        qabs = (absEff_[end]-absEff_[end-1])/(diaVal_[end]-diaVal_[end-1])*(pDia - diaVal_[end])+absEff_[end];
        qsct = (sctEff_[end]-sctEff_[end-1])/(diaVal_[end]-diaVal_[end-1])*(pDia - diaVal_[end])+sctEff_[end];
        asyF = (asyFac_[end]-asyFac_[end-1])/(diaVal_[end]-diaVal_[end-1])*(pDia - diaVal_[end])+asyFac_[end];
    }
    else
    {
        for(label i=0;i<diaVal_.size()-1;i++)
        {
            label i1 = i;
            label i2 = i+1;
            if(pDia > diaVal_[i1] && pDia <= diaVal_[i2])
            {
                scalar weight = (pDia-diaVal_[i1])/(diaVal_[i2]-diaVal_[i1]);
                qabs = (1.0-weight)*absEff_[i1]+weight*absEff_[i2];
                qsct = (1.0-weight)*sctEff_[i1]+weight*sctEff_[i2];
                asyF = (1.0-weight)*asyFac_[i1]+weight*asyFac_[i2];
            }
        }
    }

    //Info << "qabs: " << qabs << endl;
    //Info << "qsct: " << qsct << endl;

    if (nBands_ !=1)
    {
     QAbs = qabs;
     QSct = qsct; 
     asyFac = asyF;
    //Info << "qabs: " << qabs << endl;
    //Info << "qsct: " << qsct << endl;
    }
    else
    {
     List<scalar> tau(numPropBands_,0.);
     scalar tauGray = 0.;
     scalar rSum = 0.;

     scalar asyNum = 0.;
     scalar asyDen = 0.;
     
     forAll(tau,pBandI)
     { 
        tau[pBandI] = 1 - exp(-0.25*3.1415*pDia*pDia*np*(qabs[pBandI]+qsct[pBandI])*beamLen_);
        tauGray += energyFrac_[pBandI]*tau[pBandI];
        rSum += energyFrac_[pBandI]*tau[pBandI]*qabs[pBandI]/(qabs[pBandI]+qsct[pBandI]+VSMALL);
        //Info << " tau tauGray   rSum pBandI: "  << tau[pBandI]  << " " << tauGray << " " << rSum << endl;
        asyNum += energyFrac_[pBandI]*qsct[pBandI]*asyF[pBandI];
        asyDen += energyFrac_[pBandI]*qsct[pBandI];
     }
     
     scalar kextGray = -log(1. - tauGray)/(beamLen_+VSMALL);
     //Info << "kextgray: " << kextGray << endl;
     QAbs[0] = kextGray*rSum/(tauGray+VSMALL);
     QSct[0] = kextGray - QAbs[0]; 
     scalar kOverQ = 0.25*3.1415*pDia*pDia*np;
     QAbs[0] = QAbs[0]/(kOverQ + VSMALL);
     QSct[0] = QSct[0]/(kOverQ + VSMALL);

     asyFac[0] = asyNum/(asyDen + VSMALL);
    }   
    //Info << "QAbs: "  << QAbs << endl;
    //Info << "np: "   << np << endl;
    //Info << "pDia: " << pDia << endl;
    //Info << "beamLength: " << beamLen_ << endl;
    return;
}

template<class CloudType>
inline Foam::scalar Foam::ThermoCloud<CloudType>::Tmax() const
{
    scalar T = -GREAT;
    scalar n = 0;
    forAllConstIter(typename ThermoCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        T = max(T, p.T());
        n++;
    }

    reduce(T, maxOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
inline Foam::scalar Foam::ThermoCloud<CloudType>::Tmin() const
{
    scalar T = GREAT;
    scalar n = 0;
    forAllConstIter(typename ThermoCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        T = min(T, p.T());
        n++;
    }

    reduce(T, minOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setCScatCoeffs()
{
     const Foam::radiation::radiationModel& radModel_ = this->mesh().template lookupObject<Foam::radiation::radiationModel>("radiationProperties");
     const Foam::radiation::fvDOM& fvdomModel_=refCast<const Foam::radiation::fvDOM>(radModel_);

     label nRays = fvdomModel_.nRay();

     CScatCoeffs_ = 0.;
     const Foam::vector dRay0 = fvdomModel_.IRay(0).d();
     scalar omega0 = fvdomModel_.IRay(0).omega();

     for (int i=1; i < nRays; i++)
     {
       scalar s0Dotsi = dRay0 & fvdomModel_.IRay(i).d();
       scalar omegaI  = fvdomModel_.IRay(i).omega();
      
       CScatCoeffs_[0] += omegaI*s0Dotsi;
       CScatCoeffs_[1] += omegaI*s0Dotsi*s0Dotsi;
       CScatCoeffs_[2] += omegaI*s0Dotsi*s0Dotsi*s0Dotsi;
       CScatCoeffs_[3] += omegaI*s0Dotsi*s0Dotsi*s0Dotsi*s0Dotsi;
     } 

     CScatCoeffs_[0] *= 1./(4.*constant::mathematical::pi - omega0);
     CScatCoeffs_[1] *= 1./(4.*constant::mathematical::pi - omega0);
     CScatCoeffs_[2] *= 1./(4.*constant::mathematical::pi - omega0);
     CScatCoeffs_[3] *= 1./(4.*constant::mathematical::pi - omega0);
}

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setCScat()
{
     const scalar dt = this->db().time().deltaTValue();
     //Info << "Here setting cscat" << endl;
     forAll(CScat_, bandI)
     {
            scalarField& cscat = CScat_[bandI].field();
            scalarField asyFacVal = radAreaPScAsy_[bandI].field()/(radAreaPSc_[bandI].field() + VSMALL);    
            scalarField A_asy = 2.*asyFacVal/(1. + asyFacVal*asyFacVal);
            scalarField B_asy = (1. - asyFacVal*asyFacVal)/(1. + asyFacVal*asyFacVal);
            cscat =  1./( B_asy*(1. + A_asy*CScatCoeffs_[0] + A_asy*A_asy*CScatCoeffs_[1] + Foam::pow(A_asy,3.)*CScatCoeffs_[2] + Foam::pow(A_asy,4.)*CScatCoeffs_[3]) + VSMALL);
     }
}
// ************************************************************************* //

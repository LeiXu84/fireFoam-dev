/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "ParticleCollectorFmg.H"
#include "Pstream.H"
#include "surfaceWriter.H"
#include "unitConversion.H"
#include "Random.H"
#include "triangle.H"
#include "cloud.H"
#include "IOmanip.H"
#include "ListOps.H"
#include "Time.H"
#include "fvMesh.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleCollectorFmg<CloudType>::makeLogFile
(
    const faceList& faces,
    const Field<point>& points,
    const Field<scalar>& area
)
{
    // Create the output file if not already created
    if (log_)
    {
        if (debug)
        {
            Info<< "Creating output file" << endl;
        }

        if (Pstream::master())
        {
            // Create directory if does not exist
            mkDir(this->writeTimeDir());

            // Open new file at start up
            outputFilePtr_.reset
            (
                new OFstream(this->writeTimeDir()/(type() + ".dat"))
            );

            outputFilePtr_()
                << "# Type       : " << type() << nl
                << "# Source     : " << this->modelName() << nl
                << "# Bins       : " << faces.size() << nl
                << "# Total area : " << sum(area) << nl;

            switch (mode_)
            {
                case mtConcentricCircle:
                {
                    outputFilePtr_()
                        << "# nSectors   : " << nSector_ << nl
                        << "# nRadii     : " << radius_.size() << nl
                        << "# center     : " << points_[0] << nl;

                    forAll(radius_,i)
                    {
                        outputFilePtr_()
                            << "# radius " << i << ":" << radius_[i] << nl;
                    }

                    break;
                }
                default:
                {
                }
            }

            outputFilePtr_()
                << "# Geometry   :" << nl
                << '#'
                << tab << "Bin"
                << tab << "(Centre_x Centre_y Centre_z)"
                << tab << "Area"
                << nl;

            forAll(faces, i)
            {
                outputFilePtr_()
                    << '#'
                    << tab << i
                    << tab << faces[i].centre(points)
                    << tab << area[i]
                    << nl;
            }

            outputFilePtr_()
                << '#' << nl
                << "# Output format:" << nl;

            forAll(faces, i)
            {
                word id = Foam::name(i);
                word binId = "bin_" + id;

                outputFilePtr_()
                    << '#'
                    << tab << "Time"
                    << tab << binId
                    << tab << "mass[" << id << "]"
                    << tab << "massFlowRate[" << id << "]"
                    << tab << "mom[" << id << "]"
                    << tab << "momRate[" << id << "]"
                    << endl;
            }
        }
    }
}


template<class CloudType>
void Foam::ParticleCollectorFmg<CloudType>::initPolygons
(
    const List<Field<point>>& polygons
)
{
    mode_ = mtPolygon;

    label nPoints = 0;
    forAll(polygons, polyI)
    {
        label np = polygons[polyI].size();
        if (np < 3)
        {
            FatalIOErrorInFunction(this->coeffDict())
                << "polygons must consist of at least 3 points"
                << exit(FatalIOError);
        }

        nPoints += np;
    }

    label pointOffset = 0;
    points_.setSize(nPoints);
    faces_.setSize(polygons.size());
    faceTris_.setSize(polygons.size());
    area_.setSize(polygons.size());
    forAll(faces_, facei)
    {
        const Field<point>& polyPoints = polygons[facei];
        face f(identity(polyPoints.size()) + pointOffset);
        UIndirectList<point>(points_, f) = polyPoints;
        area_[facei] = f.mag(points_);

        DynamicList<face> tris;
        f.triangles(points_, tris);
        faceTris_[facei].transfer(tris);

        faces_[facei].transfer(f);

        pointOffset += polyPoints.size();
    }
}


template<class CloudType>
void Foam::ParticleCollectorFmg<CloudType>::initConcentricCircles()
{
    mode_ = mtConcentricCircle;

    vector origin(this->coeffDict().lookup("origin"));

    this->coeffDict().lookup("radius") >> radius_;
    nSector_ = readLabel(this->coeffDict().lookup("nSector"));

    label nS = nSector_;

    vector refDir;
    if (nSector_ > 1)
    {
        refDir = this->coeffDict().lookup("refDir");
        refDir -= normal_[0]*(normal_[0] & refDir);
        refDir /= mag(refDir);
    }
    else
    {
        // set 4 quadrants for single sector cases
        nS = 4;

        vector tangent = Zero;
        scalar magTangent = 0.0;

        Random rnd(1234);
        while (magTangent < SMALL)
        {
            vector v = rnd.vector01();

            tangent = v - (v & normal_[0])*normal_[0];
            magTangent = mag(tangent);
        }

        refDir = tangent/magTangent;
    }

    scalar dTheta = 5.0;
    scalar dThetaSector = 360.0/scalar(nS);
    label intervalPerSector = max(1, ceil(dThetaSector/dTheta));
    dTheta = dThetaSector/scalar(intervalPerSector);

    label nPointPerSector = intervalPerSector + 1;

    label nPointPerRadius = nS*(nPointPerSector - 1);
    label nPoint = radius_.size()*nPointPerRadius;
    label nFace = radius_.size()*nS;

    // add origin
    nPoint++;

    points_.setSize(nPoint);
    faces_.setSize(nFace);
    area_.setSize(nFace);

    coordSys_ = cylindricalCS("coordSys", origin, normal_[0], refDir, false);

    List<label> ptIDs(identity(nPointPerRadius));

    points_[0] = origin;

    // points
    forAll(radius_, radI)
    {
        label pointOffset = radI*nPointPerRadius + 1;

        for (label i = 0; i < nPointPerRadius; i++)
        {
            label pI = i + pointOffset;
            point pCyl(radius_[radI], degToRad(i*dTheta), 0.0);
            points_[pI] = coordSys_.globalPosition(pCyl);
        }
    }

    // faces
    DynamicList<label> facePts(2*nPointPerSector);
    forAll(radius_, radI)
    {
        if (radI == 0)
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                // append origin point
                facePts.append(0);

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = ptIDs.fcIndex(i - 1) + 1;
                    facePts.append(id);
                }

                label facei = secI + radI*nS;

                faces_[facei] = face(facePts);
                area_[facei] = faces_[facei].mag(points_);
            }
        }
        else
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                label offset = (radI - 1)*nPointPerRadius + 1;

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }
                for (label ptI = nPointPerSector-1; ptI >= 0; ptI--)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + nPointPerRadius + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }

                label facei = secI + radI*nS;

                faces_[facei] = face(facePts);
                area_[facei] = faces_[facei].mag(points_);
            }
        }
    }
}


template<class CloudType>
void Foam::ParticleCollectorFmg<CloudType>::collectParcelPolygon
(
    const point& p1,
    const point& p2
) const
{
    label dummyNearType = -1;
    label dummyNearLabel = -1;

    forAll(faces_, facei)
    {
        const label facePoint0 = faces_[facei][0];

        const point& pf = points_[facePoint0];

        const scalar d1 = normal_[facei] & (p1 - pf);
        const scalar d2 = normal_[facei] & (p2 - pf);

        if (sign(d1) == sign(d2))
        {
            // did not cross polygon plane
            continue;
        }

        // intersection point
        const point pIntersect = p1 + (d1/(d1 - d2))*(p2 - p1);

        const List<face>& tris = faceTris_[facei];

        // identify if point is within poly bounds
        forAll(tris, triI)
        {
            const face& tri = tris[triI];
            triPointRef t
            (
                points_[tri[0]],
                points_[tri[1]],
                points_[tri[2]]
            );

            if (t.classify(pIntersect, dummyNearType, dummyNearLabel))
            {
                hitFaceIDs_.append(facei);
            }
        }
    }
}


template<class CloudType>
void Foam::ParticleCollectorFmg<CloudType>::collectParcelConcentricCircles
(
    const point& p1,
    const point& p2
) const
{
    label secI = -1;

    const scalar d1 = normal_[0] & (p1 - coordSys_.origin());
    const scalar d2 = normal_[0] & (p2 - coordSys_.origin());

    if (sign(d1) == sign(d2))
    {
        // did not cross plane
        return;
    }

    // intersection point in cylindrical co-ordinate system
    const point pCyl = coordSys_.localPosition(p1 + (d1/(d1 - d2))*(p2 - p1));

    scalar r = pCyl[0];

    if (r < radius_.last())
    {
        label radI = 0;
        while (r > radius_[radI])
        {
            radI++;
        }

        if (nSector_ == 1)
        {
            secI = 4*radI;
        }
        else
        {
            scalar theta = pCyl[1] + constant::mathematical::pi;

            secI =
                nSector_*radI
              + floor
                (
                    scalar(nSector_)*theta/constant::mathematical::twoPi
                );
        }
    }

    if (secI != -1)
    {
        hitFaceIDs_.append(secI);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleCollectorFmg<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const Time& time = mesh.time();
    scalar timeNew = time.value();
    scalar timeElapsed = timeNew - timeOld_;

    totalTime_ += timeElapsed;

    const scalar alpha = (totalTime_ - timeElapsed)/totalTime_;
    const scalar beta = timeElapsed/(totalTime_+VSMALL);

    forAll(faces_, facei)
    {
        massTotal_[facei] += mass_[facei];
        // kvm, massFlowRate_ calculation is wrong when resetOnWrite set to false
        massFlowRate_[facei] = alpha*massFlowRate_[facei] + beta*mass_[facei]/(timeElapsed+VSMALL);
        massFlux_[facei] = massFlowRate_[facei]/area_[facei];
        momTotal_[facei] += mom_[facei];
        momRate_[facei] = alpha*momRate_[facei] + beta*mom_[facei]/(timeElapsed+VSMALL);
    }

    const label proci = Pstream::myProcNo();

    Info<< type() << " output:" << nl;

    Field<scalar> faceMassTotal(mass_.size(), 0.0);
    this->getModelProperty("massTotal", faceMassTotal);

    Field<scalar> faceMomTotal(mom_.size(), 0.0);
    this->getModelProperty("momTotal", faceMomTotal);

    Field<scalar> faceMassFlowRate(massFlowRate_.size(), 0.0);
    this->getModelProperty("massFlowRate", faceMassFlowRate);

    Field<scalar> faceMassFlux(massFlux_.size(), 0.0);
    this->getModelProperty("massFlux", faceMassFlux);

    Field<scalar> faceMomRate(momRate_.size(), 0.0);
    this->getModelProperty("momRate", faceMomRate);


    scalar sumTotalMass = 0.0;
    scalar sumTotalMom = 0.0;
    scalar sumAverageMFR = 0.0;
    scalar sumAverageMomRate = 0.0;
    forAll(faces_, facei)
    {
        scalarList allProcMass(Pstream::nProcs());
        allProcMass[proci] = massTotal_[facei];
        Pstream::gatherList(allProcMass);
        faceMassTotal[facei] += sum(allProcMass);

        scalarList allProcMom(Pstream::nProcs());
        allProcMom[proci] = momTotal_[facei];
        Pstream::gatherList(allProcMom);
        faceMomTotal[facei] += sum(allProcMom);

        scalarList allProcMassFlowRate(Pstream::nProcs());
        allProcMassFlowRate[proci] = massFlowRate_[facei];
        Pstream::gatherList(allProcMassFlowRate);
        faceMassFlowRate[facei] += sum(allProcMassFlowRate);

        scalarList allProcMassFlux(Pstream::nProcs());
        allProcMassFlux[proci] = massFlux_[facei];
        Pstream::gatherList(allProcMassFlux);
        faceMassFlux[facei] += sum(allProcMassFlux);

        scalarList allProcMomRate(Pstream::nProcs());
        allProcMomRate[proci] = momRate_[facei];
        Pstream::gatherList(allProcMomRate);
        faceMomRate[facei] += sum(allProcMomRate);

        sumTotalMass += faceMassTotal[facei];
        sumTotalMom += faceMomTotal[facei];
        sumAverageMFR += faceMassFlowRate[facei];
        sumAverageMomRate += faceMomRate[facei];

        if (outputFilePtr_.valid())
        {
            outputFilePtr_()
                << time.timeName()
                << tab << facei
                << tab << faceMassTotal[facei]
                << tab << faceMassFlowRate[facei]
                << tab << faceMomTotal[facei]
                << tab << faceMomRate[facei]
                << endl;
        }

        if (sampleParticles_)
        {
            // gather particle statistics onto master

            Pstream::gatherList(diameters_[facei]);
            Pstream::gatherList(velocityMagnitudes_[facei]);
            Pstream::gatherList(numberParticles_[facei]);
            Pstream::gatherList(positions_[facei]);

            if (Pstream::master())
            {
                // create a single list of each property
                scalarList localDiameters;
                scalarList localVelocityMagnitudes;
                scalarList localNumberParticles;
                vectorList localPositions;
                for(label proc=0;proc<Pstream::nProcs();proc++)
                {
                    localDiameters.append(diameters_[facei][proc]); 
                    localVelocityMagnitudes.append(velocityMagnitudes_[facei][proc]); 
                    localNumberParticles.append(numberParticles_[facei][proc]); 
                    localPositions.append(positions_[facei][proc]); 
                }


                fileName  basedir(this->outputDir()/time.timeName());
                fileName  filename(
                             basedir/"parcels_"
                            +Foam::name(facei)+
                            ".dat"
                          );

                mkDir(basedir);

                OFstream os(filename);

                if (!os.good())
                {
                    FatalErrorIn
                        (
                         "ParticleCollectorFmg::write()"
                        )
                        << "Cannot open file for writing " 
                        << filename
                        << exit(FatalError);
                }

                
                os << "# center " << faces_[0].centre(points_) << nl;
                os << "# d(m) u(m/s) nP x y z\n";

                forAll(localDiameters,d)
                {
                    os
                        << localDiameters[d]
                        << " "
                        << localVelocityMagnitudes[d]
                        << " "
                        << localNumberParticles[d]
                        << " "
                        << localPositions[d][0]
                        << " "
                        << localPositions[d][1]
                        << " "
                        << localPositions[d][2]
                        << nl;
                }
                os << endl;

            }
            // reset particle statistics
            diameters_[facei][Pstream::myProcNo()].clear();
            velocityMagnitudes_[facei][Pstream::myProcNo()].clear();
            numberParticles_[facei][Pstream::myProcNo()].clear();
            positions_[facei][Pstream::myProcNo()].clear();
        }
    }

    Info<< "    sum(total mass) = " << sumTotalMass << nl
        << "    sum(total momentum) = " << sumTotalMom << nl
        << "    sum(average mass flow rate) = " << sumAverageMFR << nl
        << endl;


    if (surfaceFormat_ != "none")
    {
        if (Pstream::master())
        {
            autoPtr<surfaceWriter> writer(surfaceWriter::New(surfaceFormat_));

            writer->write
            (
                this->writeTimeDir(),
                "collector",
                points_,
                faces_,
                "massTotal",
                faceMassTotal,
                false
            );

            writer->write
            (
                this->writeTimeDir(),
                "collector",
                points_,
                faces_,
                "massFlowRate",
                faceMassFlowRate,
                false
            );

            writer->write
            (
                this->writeTimeDir(),
                "collector",
                points_,
                faces_,
                "massFlux",
                faceMassFlux,
                false
            );

            writer->write
            (
                this->writeTimeDir(),
                "collector",
                points_,
                faces_,
                "momTotal",
                faceMomTotal,
                false
            );

            writer->write
            (
                this->writeTimeDir(),
                "collector",
                points_,
                faces_,
                "momRate",
                faceMomRate,
                false
            );
        }
    }


    if (resetOnWrite_)
    {
        Field<scalar> dummy(faceMassTotal.size(), 0.0);
        Field<scalar> dummy2(faceMomTotal.size(), 0.0);
        this->setModelProperty("massTotal", dummy);
        this->setModelProperty("momTotal", dummy2);
        this->setModelProperty("massFlowRate", dummy);
        this->setModelProperty("massFlux", dummy);
        this->setModelProperty("momRate", dummy);

        timeOld_ = timeNew;
        totalTime_ = 0.0;
    }
    else
    {
        this->setModelProperty("massTotal", faceMassTotal);
        this->setModelProperty("momTotal", faceMomTotal);
        this->setModelProperty("massFlowRate", faceMassFlowRate);
        this->setModelProperty("massFlux", faceMassFlux);
        this->setModelProperty("momRate", faceMomRate);
    }

    forAll(faces_, facei)
    {
        mass_[facei] = 0.0;
        mom_[facei] = 0.0;
        massTotal_[facei] = 0.0;
        momTotal_[facei] = 0.0;
        massFlowRate_[facei] = 0.0;
        massFlux_[facei] = 0.0;
        momRate_[facei] = 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleCollectorFmg<CloudType>::ParticleCollectorFmg
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    mode_(mtUnknown),
    parcelType_(this->coeffDict().lookupOrDefault("parcelType", -1)),
    removeCollected_(this->coeffDict().lookup("removeCollected")),
    points_(),
    faces_(),
    faceTris_(),
    nSector_(0),
    radius_(),
    coordSys_(false),
    normal_(),
    negateParcelsOppositeNormal_
    (
        readBool(this->coeffDict().lookup("negateParcelsOppositeNormal"))
    ),
    surfaceFormat_(this->coeffDict().lookup("surfaceFormat")),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    totalTime_(0.0),
    mass_(),
    mom_(),
    massTotal_(),
    momTotal_(),
    massFlowRate_(),
    massFlux_(),
    momRate_(),
    diameters_(),
    velocityMagnitudes_(),
    numberParticles_(),
    positions_(),
    sampleParticles_(this->coeffDict().lookupOrDefault("sampleParticles",false)),
    log_(this->coeffDict().lookup("log")),
    outputFilePtr_(),
    timeOld_(owner.mesh().time().value()),
    hitFaceIDs_()
{
    normal_ /= mag(normal_);

    word mode(this->coeffDict().lookup("mode"));
    if (mode == "polygon")
    {
        List<Field<point> > polygons(this->coeffDict().lookup("polygons"));

        initPolygons(polygons);

        vector n0(this->coeffDict().lookup("normal"));
        normal_ = vectorField(faces_.size(), n0);
    }
    else if (mode == "polygonWithNormal")
    {
        List<Tuple2<Field<point>, vector> > polygonAndNormal
        (
            this->coeffDict().lookup("polygons")
        );

        List<Field<point> > polygons(polygonAndNormal.size());
        normal_.setSize(polygonAndNormal.size());

        forAll(polygons, polyI)
        {
            polygons[polyI] = polygonAndNormal[polyI].first();
            normal_[polyI] = polygonAndNormal[polyI].second();
            normal_[polyI] /= mag(normal_[polyI]) + ROOTVSMALL;
        }

        initPolygons(polygons);
    }
    else if (mode == "concentricCircle")
    {
        vector n0(this->coeffDict().lookup("normal"));
        normal_ = vectorField(1, n0);

        initConcentricCircles();
    }
    else
    {
        FatalIOErrorIn
        (
            "Foam::ParticleCollectorFmg<CloudType>::ParticleCollectorFmg"
            "("
                "const dictionary&,"
                "CloudType&, "
                "const word&"
            ")",
            this->coeffDict()
        )
            << "Unknown mode " << mode << ".  Available options are "
            << "polygon, polygonWithNormal and concentricCircle"
            << exit(FatalIOError);
    }

    mass_.setSize(faces_.size(), 0.0);
    mom_.setSize(faces_.size(), 0.0);
    massTotal_.setSize(faces_.size(), 0.0);
    momTotal_.setSize(faces_.size(), 0.0);
    massFlowRate_.setSize(faces_.size(), 0.0);
    massFlux_.setSize(faces_.size(), 0.0);
    momRate_.setSize(faces_.size(), 0.0);

    makeLogFile(faces_, points_, area_);


    diameters_.setSize(faces_.size());
    velocityMagnitudes_.setSize(faces_.size());
    numberParticles_.setSize(faces_.size());
    positions_.setSize(faces_.size());

    forAll(faces_,facei)
    {
        diameters_[facei].setSize(Pstream::nProcs());
        velocityMagnitudes_[facei].setSize(Pstream::nProcs());
        numberParticles_[facei].setSize(Pstream::nProcs());
        positions_[facei].setSize(Pstream::nProcs());
    }

}


template<class CloudType>
Foam::ParticleCollectorFmg<CloudType>::ParticleCollectorFmg
(
    const ParticleCollectorFmg<CloudType>& pc
)
:
    CloudFunctionObject<CloudType>(pc),
    mode_(pc.mode_),
    parcelType_(pc.parcelType_),
    removeCollected_(pc.removeCollected_),
    points_(pc.points_),
    faces_(pc.faces_),
    faceTris_(pc.faceTris_),
    nSector_(pc.nSector_),
    radius_(pc.radius_),
    coordSys_(pc.coordSys_),
    normal_(pc.normal_),
    negateParcelsOppositeNormal_(pc.negateParcelsOppositeNormal_),
    surfaceFormat_(pc.surfaceFormat_),
    resetOnWrite_(pc.resetOnWrite_),
    totalTime_(pc.totalTime_),
    mass_(pc.mass_),
    mom_(pc.mom_),
    massTotal_(pc.massTotal_),
    momTotal_(pc.momTotal_),
    massFlowRate_(pc.massFlowRate_),
    massFlux_(pc.massFlux_),
    momRate_(pc.momRate_),
    diameters_(pc.diameters_),
    velocityMagnitudes_(pc.velocityMagnitudes_),
    numberParticles_(pc.numberParticles_),
    positions_(pc.positions_),
    sampleParticles_(pc.sampleParticles_),
    log_(pc.log_),
    outputFilePtr_(),
    timeOld_(0.0),
    hitFaceIDs_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleCollectorFmg<CloudType>::~ParticleCollectorFmg()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleCollectorFmg<CloudType>::postMove
(
    parcelType& p,
    const label celli,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    if ((parcelType_ != -1) && (parcelType_ != p.typeId()))
    {
        return;
    }

    // slightly extend end position to avoid falling within tracking tolerances
    const point position1 = position0 + 1.0001*(p.position() - position0);

    hitFaceIDs_.clear();

    switch (mode_)
    {
        case mtPolygon:
        {
            collectParcelPolygon(position0, position1);
            break;
        }
        case mtConcentricCircle:
        {
            collectParcelConcentricCircles(position0, position1);
            break;
        }
        default:
        {
        }
    }


    forAll(hitFaceIDs_, i)
    {
        label facei = hitFaceIDs_[i];
        scalar m = p.nParticle()*p.mass();
        scalar d = p.d();
        scalar u = mag(p.U());
        /*scalar Uy = p.U().component(1);*/
        /*scalar Uz = p.U().component(2);*/
        scalar np = p.nParticle();
        vector position = p.position();

        /*cout << "diameter " << d << " velocity " << u << " nParticle " << np << "\n";*/

        if(sampleParticles_)
        {
            diameters_[facei][Pstream::myProcNo()].append(d);
            velocityMagnitudes_[facei][Pstream::myProcNo()].append(u);
            numberParticles_[facei][Pstream::myProcNo()].append(np);
            positions_[facei][Pstream::myProcNo()].append(position);
        }

        if (negateParcelsOppositeNormal_)
        {
            vector Uhat = p.U();
            Uhat /= mag(Uhat) + ROOTVSMALL;
            if ((Uhat & normal_[facei]) < 0)
            {
                m *= -1.0;
            }
        }

        // add mass contribution
        mass_[facei] += m; // should be *0.25 if nSector_ == 1
        vector Uhat = p.U();
        scalar Unormal = Uhat & normal_[facei];
        mom_[facei] += m*Unormal; // should be *0.25 if nSector_ == 1

        if (nSector_ == 1)
        {
            mass_[facei + 1] += m; // should be *0.25
            mass_[facei + 2] += m; // should be *0.25
            mass_[facei + 3] += m; // should be *0.25
            mom_[facei + 1] += m*Unormal; // should be *0.25
            mom_[facei + 2] += m*Unormal; // should be *0.25
            mom_[facei + 3] += m*Unormal; // should be *0.25
        }

        if (removeCollected_)
        {
            keepParticle = false;
        }
    }
}


// ************************************************************************* //

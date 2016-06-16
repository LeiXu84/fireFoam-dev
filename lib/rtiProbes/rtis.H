/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::rtis

Group
    grpUtilitiesFunctionObjects

Description
    Set of locations to sample.

    Call write() to sample and write files.

SourceFiles
    rtis.C

\*---------------------------------------------------------------------------*/

#ifndef rtis_H
#define rtis_H

#include "functionObject.H"
#include "Switch.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "pointField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceMesh.H"
#include "wordReList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DEBUG(x) std::cout.precision(10); std::cout.setf(std::ios::fixed,std::ios::floatfield);std::cout << "["<< __FILE__ << ":" << __LINE__ << "] p" << Pstream::myProcNo() << " " << #x " = " << x << std::endl;
#define TRACE(s) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] p" << Pstream::myProcNo() << " " << #s << std::endl; s;

namespace Foam
{

// Forward declaration of classes
class Time;
class objectRegistry;
class dictionary;
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class rtis Declaration
\*---------------------------------------------------------------------------*/

class rtis
:
    public functionObject,
    public pointField
{
protected:

    // Protected classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}
        };


    // Private data

        //- Name of this set of rtis,
        //  Also used as the name of the rtis directory.
        word name_;

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;


        // Read from dictonary

            //- Names of fields to rti
            wordReList fieldSelection_;

            //- Fixed locations, default = yes
            //  Note: set to false for moving mesh calations where locations
            //        should move with the mesh
            bool fixedLocations_;

            //- Interpolation scheme name
            //  Note: only possible when fixedLocations_ is true
            word interpolationScheme_;

            //- Flag for determining restart condition
            Switch cleanRestart_;


        // Calculated

            //- Categorized scalar/vector/tensor vol fields
            fieldGroup<scalar> scalarFields_;
            fieldGroup<vector> vectorFields_;
            fieldGroup<sphericalTensor> sphericalTensorFields_;
            fieldGroup<symmTensor> symmTensorFields_;
            fieldGroup<tensor> tensorFields_;

            //- Categorized scalar/vector/tensor surf fields
            fieldGroup<scalar> surfaceScalarFields_;
            fieldGroup<vector> surfaceVectorFields_;
            fieldGroup<sphericalTensor> surfaceSphericalTensorFields_;
            fieldGroup<symmTensor> surfaceSymmTensorFields_;
            fieldGroup<tensor> surfaceTensorFields_;

            // Cells to be rtid (obtained from the locations)
            labelList elementList_;

            // Faces to be rtid
            labelList faceList_;

            //- Current open files
            HashPtrTable<OFstream> rtiFilePtrs_;

            //-RTI related
            scalar RTI_;
            scalar C_;
            scalarList initialTemperatures_;
            scalar activationTemperature_;
            scalarList rtiValues_;
            List < scalar > dTeOld_;

    // Private Member Functions

        //- Clear old field groups
        void clearFieldGroups();

        //- Append fieldName to the appropriate group
        label appendFieldGroup(const word& fieldName, const word& fieldType);

        //- Classify field types, returns the number of fields
        label classifyFields();

        //- Find cells and faces containing rtis
        virtual void findElements(const fvMesh&);

        //- Classify field type and Open/close file streams,
        //  returns number of fields to sample
        label prepare();

        void prepareRTI();

        void writeRTI();

        //- Compute link temperature based on RTI
        scalar calculateLinkTemperature(label rtiI, scalar T, scalar U);

private:

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        //- Sample and write all the fields of the given type
        template<class Type>
        void sampleAndWrite(const fieldGroup<Type>&);

        //- Sample and write all the surface fields of the given type
        template<class Type>
        void sampleAndWriteSurfaceFields(const fieldGroup<Type>&);

        //- Disallow default bitwise copy construct
        rtis(const rtis&);

        //- Disallow default bitwise assignment
        void operator=(const rtis&);


public:

    //- Runtime type information
    TypeName("rtis");


    // Constructors

        //- Construct from Time and dictionary
        rtis
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        rtis
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~rtis();


    // Member Functions

        //- Return names of fields to rti
        virtual const wordReList& fieldNames() const
        {
            return fieldSelection_;
        }

        //- Return locations to rti
        virtual const pointField& rtiLocations() const
        {
            return *this;
        }

        //- Return location for rti i
        virtual const point& rti(const label i) const
        {
            return operator[](i);
        }

        //- Cells to be rtid (obtained from the locations)
        const labelList& elements() const
        {
            return elementList_;
        }

        //- Read the rtis
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet();

        //- Sample and write
        virtual bool write();

        //- Read probe file for restart
        typedef Tuple2<scalarList,DynamicList<scalarList> > probeData;
        
        probeData readProbeFile(std::string testFile);

        bool interpolateProbeData(const probeData& data);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state)
        {}

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Sample a single vol field on all sample locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Sample a single scalar field on all sample locations
        template<class Type>
        tmp<Field<Type>> sampleSurfaceFields(const word& fieldName) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "rtisTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM and modified by the
    authors described below.

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
    Foam::searchableCone

Description
    Surface geometry with a cone shape, which can be used with
    snappyHexMesh.

Authors
    Stanislau Stasheuski, Aalto University, 2023
    stanislau.stasheuski@aalto.fi

    Mark Olesen, ESI, 2020.

Usage
    \table
        Property  | Description                      | Required | Default
        point1    | Point at one end of cone axis    | yes      |
        point2    | Point at other end of cone axis  | yes      |
        radius1   | Cone radius at one end           | yes      |
        radius2   | Cone radius at other end         | no       | 0
    \endtable

    Example specification in snappyHexMeshDict/geometry:
    \verbatim
        type       searchableCone;
        point1     (0 0 0);
        point2     (0 0 1);
        radius1    0.25;
    \endverbatim

SourceFiles
    searchableCone.C

\*---------------------------------------------------------------------------*/

#ifndef searchableCone_H
#define searchableCone_H

#include "treeBoundBox.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                     Class searchableCone Declaration
\*---------------------------------------------------------------------------*/

class searchableCone
:
    public searchableSurface
{
    // Private Member Data

        //- 'left' point
        const point point1_;

        //- 'right' point
        const point point2_;

        //- Length of vector point2-point1
        const scalar magDir_;

        //- Normalised vector point2-point1
        const vector unitDir_;

        //- Radius squared at point1_
        const scalar radius1_;

        //- Radius squared at point2_
        const scalar radius2_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        //- Find nearest point and normal on cone.
        void findNearestAndNormal
        (
            const point& sample,
            const scalar nearestDistSqr,
            pointIndexHit& info,
            vector& nearNormal
        ) const;


        scalar radius2(const point& pt) const;

        //- Find intersection with cone
        void findLineAll
        (
            const point& start,
            const point& end,
            pointIndexHit& near,
            pointIndexHit& far
        ) const;

        //- Return the boundBox of the cone
        boundBox calcBounds() const;


public:

    //- Runtime type information
    TypeName("searchableCone");


    // Constructors

        //- Construct from components
        searchableCone
        (
            const IOobject& io,
            const point&,
            const point&,
            const scalar radius1,
            const scalar radius2
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableCone
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchableCone(const searchableCone&) = delete;


    //- Destructor
    virtual ~searchableCone();


    // Member Functions

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return true;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const
        {
            NotImplemented;

            return false;
        }


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableCone&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

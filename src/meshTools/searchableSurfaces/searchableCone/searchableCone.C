/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM but modified by authors
    described below.

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

Authors
    Stanislau Stasheuski, Aalto University, 2023
    stanislau.stasheuski@aalto.fi

    Mark Olesen, ESI, 2020.

\*---------------------------------------------------------------------------*/


#include "searchableCone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableCone, 0);
    addToRunTimeSelectionTable(searchableSurface, searchableCone, dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::searchableCone::coordinates() const
{
    tmp<pointField> tCtrs(new pointField(1, 0.5*(point1_ + point2_)));

    return tCtrs;
}


void Foam::searchableCone::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(1);
    centres[0] = 0.5*(point1_ + point2_);

    radiusSqr.setSize(1);
    radiusSqr[0] = 
        (radius1_ > radius2_)
            ? Foam::magSqr(point1_-centres[0]) + Foam::sqr(radius1_)
            : Foam::magSqr(point2_-centres[0]) + Foam::sqr(radius2_);

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(small);
}


Foam::tmp<Foam::pointField> Foam::searchableCone::points() const
{
    tmp<pointField> tPts(new pointField(2));
    pointField& pts = tPts.ref();

    pts[0] = point1_;
    pts[1] = point2_;

    return tPts;
}


void Foam::searchableCone::findNearestAndNormal
(
    const point& sample,
    const scalar nearestDistSqr,
    pointIndexHit& info,
    vector& nearNormal
) const
{
    vector v(sample - point1_);

    // Decompose sample-point1 into normal and parallel component
    const scalar parallel = (v & unitDir_);

    // Remove the parallel component and normalise
    v -= parallel*unitDir_;

    const scalar magV = mag(v);
    v /= mag(v);

    // Nearest and normal on disk at point1
    point disk1Point(point1_ + min(magV, radius1_)*v);
    vector disk1Normal(-unitDir_);

    // Nearest and normal on disk at point2
    point disk2Point(point2_ + min(magV, radius2_)*v);
    vector disk2Normal(unitDir_);

    // Nearest and normal on cone. Initialise to far-away point so if not
    // set it picks one of the disk points
    point nearCone(point::uniform(GREAT));
    vector normalCone(1, 0, 0);

    point projPt1 = point1_ + radius1_*v;
    point projPt2 = point2_ + radius2_*v;

    vector p1 = (projPt1 - point1_);
    if (mag(p1) > rootVSmall)
    {
        p1 /= mag(p1);

        // Find vector along the two end of cone
        const vector b = normalised(projPt2 - projPt1);

        // Find the vector along sample pt and pt at one end of cone
        vector a = (sample - projPt1);

        if (mag(a) <= rootVSmall)
        {
            // Exception: sample on disk1. Redo with projPt2.
            a = (sample - projPt2);

            // Find normal unitvector
            nearCone = (a & b)*b + projPt2;

            vector b1 = (p1 & b)*b;
            normalCone = normalised(p1 - b1);
        }
        else
        {
            // Find nearest point on cone surface
            nearCone = (a & b)*b + projPt1;

            // Find projection along surface of cone
            vector b1 = (p1 & b)*b;
            normalCone = normalised(p1 - b1);
        }
    }

    // Select nearest out of the 3 points (cone, disk1, disk2)

    FixedList<scalar, 3> dist;
    dist[0] = magSqr(nearCone - sample);
    dist[1] = magSqr(disk1Point - sample);
    dist[2] = magSqr(disk2Point - sample);

    const label minI = findMin(dist);

    // Snap the point to the corresponding surface

    if (minI == 0)  // Near the cone
    {
        // Closest to (infinite) cone. See if needs clipping to end disks
        {
            vector v1(nearCone - point1_);
            scalar para = (v1 & unitDir_);
            // Remove the parallel component and normalise
            v1 -= para*unitDir_;
            const scalar magV1 = mag(v1);
            v1 = v1/magV1;

            if (para < 0.0 && magV1 >= radius1_)
            {
                // Near point 1. Set point to intersection of disk and cone.
                // Keep normal from cone.
                nearCone = disk1Point;
            }
            else if (para < 0.0 && magV1 < radius1_)
            {
                // On disk1
                nearCone = disk1Point;
                normalCone = disk1Normal;
            }
            else if (para > magDir_ && magV1 >= radius2_)
            {
                // Near point 2. Set point to intersection of disk and cone.
                // Keep normal from cone.
                nearCone = disk2Point;
            }
            else if (para > magDir_ && magV1 < radius2_)
            {
                // On disk2
                nearCone = disk2Point;
                normalCone = disk2Normal;
            }
        }
        info.setPoint(nearCone);
        nearNormal = normalCone;
    }
    else if (minI == 1) // Near to disk1
    {
        info.setPoint(disk1Point);
        nearNormal = disk1Normal;
    }
    else if (minI == 2) // Near to disk2
    {
        info.setPoint(disk2Point);
        nearNormal = disk2Normal;
    }

    if (magSqr(sample - info.rawPoint()) < nearestDistSqr)
    {
        info.setHit();
        info.setIndex(0);
    }
}


Foam::scalar Foam::searchableCone::radius2(const point& pt) const
{
    const vector x = (pt-point1_) ^ unitDir_;
    return x&x;
}


// From http://www.gamedev.net/community/forums/topic.asp?topic_id=467789 -
// intersection of cylinder with ray
void Foam::searchableCone::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    vector point1Start(start-point1_);
    vector point2Start(start-point2_);
    vector point1End(end-point1_);

    // Quick rejection of complete vector outside endcaps
    scalar s1 = point1Start&unitDir_;
    scalar s2 = point1End&unitDir_;

    if ((s1 < 0 && s2 < 0) || (s1 > magDir_ && s2 > magDir_))
    {
        return;
    }

    // Line as P = start+t*V  where V is unit vector and t=[0..mag(end-start)]
    vector V(end-start);
    scalar magV = mag(V);
    if (magV < rootVSmall)
    {
        return;
    }
    V /= magV;


    // We now get the nearest intersections to start. This can either be
    // the intersection with the end plane or with the cylinder side.

    // Get the two points (expressed in t) on the end planes. This is to
    // clip any cylinder intersection against.
    scalar tPoint1;
    scalar tPoint2;

    // Maintain the two intersections with the endcaps
    scalar tNear = vGreat;
    scalar tFar = vGreat;

    // scalar radius;
    {
        scalar s = (V&unitDir_);
        if (mag(s) > vSmall)
        {
            tPoint1 = -s1/s;
            tPoint2 = -(point2Start&unitDir_)/s;
            if (tPoint2 < tPoint1)
            {
                Swap(tPoint1, tPoint2);
            }
            if (tPoint1 > magV || tPoint2 < 0)
            {
                return;
            }

            // See if the points on the endcaps are actually inside the cylinder
            if (tPoint1 >= 0 && tPoint1 <= magV)
            {
                scalar radius = 
                (
                    mag((start+tPoint1*V-point1_)&unitDir_)
                  > mag((start+tPoint1*V-point2_)&unitDir_)
                ) 
                    ? radius2_
                    : radius1_;

                if (radius2(start+tPoint1*V) <= sqr(radius))
                {
                    tNear = tPoint1;
                }
            }
            if (tPoint2 >= 0 && tPoint2 <= magV)
            {
                scalar radius = 
                    (
                        mag((start+tPoint2*V-point1_)&unitDir_)
                      > mag((start+tPoint2*V-point2_)&unitDir_)
                    )
                        ? radius2_
                        : radius1_;

                if (radius2(start+tPoint2*V) <= sqr(radius))
                {
                    // Check if already have a near hit from point1
                    if (tNear <= 1)
                    {
                        tFar = tPoint2;
                    }
                    else
                    {
                        tNear = tPoint2;
                    }
                }
            }
        }
        else
        {
            // Vector perpendicular to cylinder. Check for outside already done
            // above so just set tpoint to allow all.
            tPoint1 = -vGreat;
            tPoint2 = vGreat;
        }
    }

    // Second order equation of the form a*t^2 + b*t + c
    scalar a, b, c;

    scalar deltaRadius = radius2_ - radius1_;
    if (mag(deltaRadius) <= rootVSmall)
    {
        vector point1Start(start - point1_);
        const vector x = point1Start ^ unitDir_;
        const vector y = V ^ unitDir_;
        const scalar d = sqr(0.5*(radius1_ + radius2_));

        a = (y&y);
        b = 2*(x&y);
        c = (x&x) - d;
    }
    else
    {
        vector v = normalised(end - start);
        scalar p  = unitDir_ & v;
        vector a1 = v - p*unitDir_;

        // Determine the end point of the cone
        point pa =
            unitDir_*radius1_*mag(point2_-point1_)/(-deltaRadius)
          + point1_;

        scalar l2 = sqr(deltaRadius) + sqr(magDir_);
        scalar sqrCosAlpha = sqr(magDir_)/l2;
        scalar sqrSinAlpha = sqr(deltaRadius)/l2;

        vector delP(start - pa);
        vector p1 = (delP - (delP&unitDir_)*unitDir_);

        a = sqrCosAlpha*((v - p*unitDir_)&(v - p*unitDir_))-sqrSinAlpha*sqr(p);
        b =
            2.0*sqrCosAlpha*(a1&p1)
          - 2.0*sqrSinAlpha*(v&unitDir_)*(delP&unitDir_);
        c =
            sqrCosAlpha
           *(
                (delP-(delP&unitDir_)*unitDir_)
              & (delP-(delP&unitDir_)*unitDir_)
            )
          - sqrSinAlpha*sqr(delP&unitDir_);
    }

    const scalar disc = sqr(b) - 4*a*c;

    scalar t1 = -vGreat;
    scalar t2 = vGreat;

    if (disc < 0)
    {
        // Fully outside
        return;
    }
    else if (disc < rootVSmall)
    {
        // Single solution
        if (mag(a) > rootVSmall)
        {
            t1 = -b/(2*a);

            if (t1 >= 0 && t1 <= magV && t1 >= tPoint1 && t1 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t1 < tNear)
                {
                    tFar = tNear;
                    tNear = t1;
                }
                else if (t1 < tFar)
                {
                    tFar = t1;
                }
            }
            else
            {
                return;
            }
        }
        else
        {
            // Aligned with axis. Check if outside radius
            if (c > 0)
            {
                return;
            }
        }
    }
    else
    {
        if (mag(a) > rootVSmall)
        {
            scalar sqrtDisc = sqrt(disc);

            t1 = (-b - sqrtDisc)/(2*a);
            t2 = (-b + sqrtDisc)/(2*a);
            if (t2 < t1)
            {
                Swap(t1, t2);
            }

            if (t1 >= 0 && t1 <= magV && t1 >= tPoint1 && t1 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t1 < tNear)
                {
                    tFar = tNear;
                    tNear = t1;
                }
                else if (t1 < tFar)
                {
                    tFar = t1;
                }
            }
            if (t2 >= 0 && t2 <= magV && t2 >= tPoint1 && t2 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t2 < tNear)
                {
                    tFar = tNear;
                    tNear = t2;
                }
                else if (t2 < tFar)
                {
                    tFar = t2;
                }
            }
        }
        else
        {
            // Aligned with axis. Check if outside radius
            if (c > 0)
            {
                return;
            }
        }
    }

    // Check tNear, tFar
    if (tNear >= 0 && tNear <= magV)
    {
        near.setPoint(start+tNear*V);
        near.setHit();
        near.setIndex(0);

        if (tFar <= magV)
        {
            far.setPoint(start+tFar*V);
            far.setHit();
            far.setIndex(0);
        }
    }
    else if (tFar >= 0 && tFar <= magV)
    {
        near.setPoint(start+tFar*V);
        near.setHit();
        near.setIndex(0);
    }
}


Foam::boundBox Foam::searchableCone::calcBounds() const
{

    // Adapted from
    // http://www.gamedev.net/community/forums
    //       /topic.asp?topic_id=338522&forum_id=20&gforum_id=0

    // Let cylinder have end points A,B and radius r,

    // Bounds in direction X (same for Y and Z) can be found as:
    // Let A.X<B.X (otherwise swap points)
    // Good approximate lowest bound is A.X-r and highest is B.X+r (precise for
    // capsule). At worst, in one direction it can be larger than needed by 2*r.

    // Accurate bounds for cylinder is
    // A.X-kx*r, B.X+kx*r
    // where
    // kx=sqrt(((A.Y-B.Y)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))

    // similar thing for Y and Z
    // (i.e.
    // ky=sqrt(((A.X-B.X)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // kz=sqrt(((A.X-B.X)^2+(A.Y-B.Y)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // )

    // How derived: geometric reasoning. Bounds of cylinder is same as for 2
    // circles centered on A and B. This sqrt thingy gives sine of angle between
    // axis and direction, used to find projection of radius.

    vector kr
    (
        sqrt(sqr(unitDir_.y()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.y()))
    );

    kr *= (radius1_ > radius2_) ? radius1_ : radius2_;

    point min = point1_ - kr;
    point max = point1_ + kr;

    min = ::Foam::min(min, point2_ - kr);
    max = ::Foam::max(max, point2_ + kr);

    return boundBox(min, max);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableCone::searchableCone
(
    const IOobject& io,
    const point& point1,
    const point& point2,
    const scalar radius1,
    const scalar radius2
)
:
    searchableSurface(io),
    point1_(point1),
    point2_(point2),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_),
    radius1_(radius1),
    radius2_(radius2)
{
    bounds() = calcBounds();
}


Foam::searchableCone::searchableCone
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableCone(
        io,
        dict.lookup("point1"),
        dict.lookup("point2"),
        dict.lookup<scalar>("radius1"),
        dict.lookupOrDefault<scalar>("radius2", small)
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableCone::~searchableCone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableCone::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchableCone::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        vector normal;
        findNearestAndNormal(samples[i], nearestDistSqr[i], info[i], normal);
    }
}


void Foam::searchableCone::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected take second one.
        pointIndexHit b;
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableCone::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        // Discard far intersection
        pointIndexHit b;
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableCone::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;
        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].setSize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].setSize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].setSize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }
}


void Foam::searchableCone::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableCone::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            pointIndexHit nearInfo;
            findNearestAndNormal
            (
                info[i].hitPoint(),
                Foam::sqr(GREAT),
                nearInfo,
                normal[i]
            );
        }
    }
}


void Foam::searchableCone::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());
    volType = volumeType::inside;

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        vector v(pt - point1_);

        // Decompose sample-point1 into normal and parallel component
        scalar parallel = v & unitDir_;

        if (parallel < 0 || parallel > magDir_)
        {
            // left or right of point1 or point2 endcaps accordingly
            volType[pointi] = volumeType::outside;
        }
        else
        {
            // Remove the parallel component
            v -= parallel*unitDir_;

            volType[pointi] = 
                (mag(v) > radius1_ + parallel * (radius2_-radius1_)/magDir_) 
                    ? volumeType::outside
                    : volumeType::inside;
        }
    }
}


// ************************************************************************* //

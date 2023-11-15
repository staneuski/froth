/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

Application
    wedgeMesh

Description
    Create a wedge from a 2D cartesian mesh.

Authors
    Stanislau Stasheuski, Aalto University, 2023
    stanislau.stasheuski@aalto.fi

    Kristjan Krebelj, University of Ljubljana, 2020

    Henry Weller, CFD-Direct, 2022.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "unitConversion.H"
#include "emptyPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label parseComponent(const word& axisName) {
    if (axisName == "X")
    {
        return vector::X;
    }
    else if (axisName == "Y")
    {
        return vector::Y;
    }
    else if (axisName == "Z")
    {
        return vector::Z;
    }

    FatalErrorInFunction
        << "Invalid axis name.\nValid names are: X, Y, Z"
        << exit(FatalError);
    return -1;  // prevent comilation warning
}

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"

    argList::addNote
    (
        "Creates a wedge from a 2D cartesian mesh by specifying normals to wedge plane\n"
        "and axis patch."
    );

    argList::validArgs.append("plane");
    argList::validArgs.append("axis normal");
    argList::validArgs.append("angle[0-5]");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedPolyMesh.H"

    const label plane = parseComponent(args.argRead<word>(1)),
                axis = parseComponent(args.argRead<word>(2));
    const scalar wedgeAngle = args.argRead<scalar>(3);

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    const fileName meshRegionSubDir = meshRegionName != polyMesh::defaultRegion
        ? meshRegionName/polyMesh::meshSubDir
        : fileName(polyMesh::meshSubDir);

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshRegionSubDir, "points"),
            meshRegionSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const point midPoint = gAverage(points);

    Info<< "Wedge: " << wedgeAngle << " deg" << nl << endl;

    const scalar halfCos = Foam::cos(degToRad(0.5*wedgeAngle)),
                 halfSin = Foam::sin(degToRad(0.5*wedgeAngle));

    forAll(points, pointI)
    {
        points[pointI].replace(
            plane,
            points[pointI].component(plane) < midPoint.component(plane)
            ? - halfSin*points[pointI].component(axis)
            : + halfSin*points[pointI].component(axis)
        );
        points[pointI].replace(
            axis,
            halfCos*points[pointI].component(axis)
        );
    }

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

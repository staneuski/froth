/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM and modified by the authors
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0;
scalar acousticCoNum = 0;

{
    const scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());

    CoNum = 0.5*gMax(sumPhi/mesh.V())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V()))*runTime.deltaTValue();

    acousticCoNum = 0.5*gMax
    (
        fvc::surfaceSum
        (
            fvc::interpolate(scalar(1)/sqrt(mixture.psi()))*mesh.magSf()
        )().primitiveField()/mesh.V()
    )*runTime.deltaTValue();

    Info<< "phi Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " acoustic max: " << acousticCoNum
        << endl;
}

// ************************************************************************* //

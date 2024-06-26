/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Class
    Foam::compressibilityModels::Chung

Authors
    Henry Weller, CFD-Direct, 2021.

Description
    Chung compressibility model.

SourceFiles
    Chung.C

\*---------------------------------------------------------------------------*/

#ifndef Chung_H
#define Chung_H

#include "barotropicCompressibilityModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressibilityModels
{

/*---------------------------------------------------------------------------*\
                            Class Chung Declaration
\*---------------------------------------------------------------------------*/

class Chung
:
    public barotropicCompressibilityModel
{
    // Private Data

        dimensionedScalar pSat_;
        dimensionedScalar psiv_;
        dimensionedScalar psil_;
        dimensionedScalar rholSat_;


public:

    //- Runtime type information
    TypeName("Chung");


    // Constructors

        //- Construct from components
        Chung
        (
            const dictionary& compressibilityProperties,
            const volScalarField& gamma,
            const word& psiName = "psi"
        );


    //- Destructor
    ~Chung()
    {}


    // Member Functions

        //- Correct the Chung compressibility
        void correct();

        //- Read phaseProperties dictionary
        bool read(const dictionary& compressibilityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressibilityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd
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

#include "NewtonianAbsorption.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(NewtonianAbsorption, 0);
    addToRunTimeSelectionTable(viscosityModel, NewtonianAbsorption, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::NewtonianAbsorption::NewtonianAbsorption
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    nu0_("nu", dimViscosity, viscosityProperties_),
    nuD_("nuD", dimViscosity, viscosityProperties_),
    x1_(viscosityProperties_.get<scalar>("x1")),
    x2_(viscosityProperties_.get<scalar>("x2")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        nu0_
    )
{
    forAll(U_.mesh().cells(),celli)
    {
        if ((U_.mesh().C()[celli][0] >= x1_) &&
            (U_.mesh().C()[celli][0] <= x2_))
        {
            nu_[celli] = nu_[celli] + nuD_.value()
                       * (U_.mesh().C()[celli][0] - x1_) / (x2_ - x1_);
        }
        else if (U_.mesh().C()[celli][0] > x2_)
        {
            nu_[celli] = nu_[celli] + nuD_.value();
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::NewtonianAbsorption::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    viscosityProperties_.readEntry("nu", nu0_);
    viscosityProperties_.readEntry("nuD", nuD_);
    viscosityProperties_.readEntry("x1", x1_);
    viscosityProperties_.readEntry("x2", x2_);
    nu_ = nu0_;

    return true;
}


// ************************************************************************* //

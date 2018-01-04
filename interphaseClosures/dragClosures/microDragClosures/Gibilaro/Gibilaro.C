/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Contributors
    Federico Municchi, TUGraz, 2017
\*---------------------------------------------------------------------------*/

#include "Gibilaro.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace microDragClosures
{
    defineTypeNameAndDebug(Gibilaro, 0);
    addToRunTimeSelectionTable(microDragClosure, Gibilaro, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::microDragClosures::Gibilaro::Gibilaro
(
    const dictionary& dict,
    const phasePair& pair
)
:
    microDragClosure(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::microDragClosures::Gibilaro::~Gibilaro()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::microDragClosures::Gibilaro::CdRe() const
{
    volScalarField alpha2
    (
        max(scalar(1) - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    return
        (4/3)
       *(17.3/alpha2 + 0.336*pair_.Re())
       *max(pair_.continuous(), pair_.continuous().residualAlpha())
       *pow(alpha2, -2.8);
}


// ************************************************************************* //

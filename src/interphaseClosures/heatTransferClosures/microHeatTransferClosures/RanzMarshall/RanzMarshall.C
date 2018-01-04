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

#include "RanzMarshall.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace microHeatTransferClosures
{
    defineTypeNameAndDebug(RanzMarshall, 0);
    addToRunTimeSelectionTable(microHeatTransferClosure, RanzMarshall, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::microHeatTransferClosures::RanzMarshall::RanzMarshall
(
    const dictionary& dict,
    const phasePair& pair
)
:
    microHeatTransferClosure(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::microHeatTransferClosures::RanzMarshall::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::microHeatTransferClosures::RanzMarshall::K() const
{
    volScalarField Nu(scalar(2) + 0.6*sqrt(pair_.Re())*cbrt(pair_.Pr()));

    return
        6.0
       *max(pair_.dispersed(), residualAlpha_)
       *pair_.continuous().kappa()
       *Nu
       /sqr(pair_.dispersed().d());
}


// ************************************************************************* //

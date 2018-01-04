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

#include "SyamlalOBrien.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace microDragClosures
{
    defineTypeNameAndDebug(SyamlalOBrien, 0);
    addToRunTimeSelectionTable(microDragClosure, SyamlalOBrien, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::microDragClosures::SyamlalOBrien::SyamlalOBrien
(
    const dictionary& dict,
    const phasePair& pair
)
:
    microDragClosure(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::microDragClosures::SyamlalOBrien::~SyamlalOBrien()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::microDragClosures::SyamlalOBrien::CdRe() const
{
    volScalarField alpha2
    (
        max(scalar(1) - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    volScalarField A(pow(alpha2, 4.14));
    volScalarField B
    (
        neg(alpha2 - 0.85)*(0.8*pow(alpha2, 1.28))
      + pos(alpha2 - 0.85)*(pow(alpha2, 2.65))
    );
    volScalarField Re(pair_.Re());
    volScalarField Vr
    (
        0.5
       *(
            A - 0.06*Re + sqrt(sqr(0.06*Re) + 0.12*Re*(2.0*B - A) + sqr(A))
        )
    );
    volScalarField CdsRe(sqr(0.63*sqrt(Re) + 4.8*sqrt(Vr)));

    return
        CdsRe
       *max(pair_.continuous(), pair_.continuous().residualAlpha())
       /sqr(Vr);
}


// ************************************************************************* //

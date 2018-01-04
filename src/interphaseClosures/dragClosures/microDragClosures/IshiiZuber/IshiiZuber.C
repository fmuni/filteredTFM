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

#include "IshiiZuber.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace microDragClosures
{
    defineTypeNameAndDebug(IshiiZuber, 0);
    addToRunTimeSelectionTable(microDragClosure, IshiiZuber, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::microDragClosures::IshiiZuber::IshiiZuber
(
    const dictionary& dict,
    const phasePair& pair
)
:
    microDragClosure(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::microDragClosures::IshiiZuber::~IshiiZuber()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::microDragClosures::IshiiZuber::CdRe() const
{
    volScalarField Re(pair_.Re());
    volScalarField Eo(pair_.Eo());

    volScalarField mud(pair_.dispersed().mu());
    volScalarField muc(pair_.continuous().mu());

    volScalarField muStar((mud + 0.4*muc)/(mud + muc));

    volScalarField muMix
    (
        muc
       *pow(max(1 - pair_.dispersed(), scalar(1e-3)), -2.5*muStar)
    );

    volScalarField ReM(Re*muc/muMix);
    volScalarField CdRe
    (
        pos(1000 - ReM)*24.0*(scalar(1) + 0.15*pow(ReM, 0.687))
      + neg(1000 - ReM)*0.44*ReM
    );

    volScalarField F((muc/muMix)*sqrt(1 - pair_.dispersed()));
    F.max(1e-3);

    volScalarField Ealpha((1 + 17.67*pow(F, 0.8571428))/(18.67*F));

    volScalarField CdReEllipse(Ealpha*0.6666*sqrt(Eo)*Re);

    return
        pos(CdReEllipse - CdRe)
       *min(CdReEllipse, Re*sqr(1 - pair_.dispersed())*2.66667)
      + neg(CdReEllipse - CdRe)*CdRe;
}


// ************************************************************************* //

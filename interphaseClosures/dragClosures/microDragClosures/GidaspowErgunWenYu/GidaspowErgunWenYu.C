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

#include "GidaspowErgunWenYu.H"
#include "phasePair.H"
#include "Ergun.H"
#include "WenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace microDragClosures
{
    defineTypeNameAndDebug(GidaspowErgunWenYu, 0);
    addToRunTimeSelectionTable(microDragClosure, GidaspowErgunWenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::microDragClosures::GidaspowErgunWenYu::GidaspowErgunWenYu
(
    const dictionary& dict,
    const phasePair& pair
)
:
    microDragClosure(dict, pair),
    Ergun_
    (
        new Ergun
        (
            dict,
            pair
        )
    ),
    WenYu_
    (
        new WenYu
        (
            dict,
            pair
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::microDragClosures::GidaspowErgunWenYu::~GidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::microDragClosures::GidaspowErgunWenYu::CdRe() const
{
    return
        pos(pair_.continuous() - 0.8)*WenYu_->CdRe()
      + neg(pair_.continuous() - 0.8)*Ergun_->CdRe();
}


// ************************************************************************* //

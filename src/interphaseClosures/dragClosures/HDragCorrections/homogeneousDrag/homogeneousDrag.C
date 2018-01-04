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

#include "homogeneousDrag.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HDragCorrections
{
    defineTypeNameAndDebug(HomogeneousDrag, 0);
    addToRunTimeSelectionTable(HDragCorrection, HomogeneousDrag, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HDragCorrections::HomogeneousDrag::HomogeneousDrag
(
    const dictionary& dict,
    const phasePair& pair
)
:
    HDragCorrection(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HDragCorrections::HomogeneousDrag::~HomogeneousDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::HDragCorrections::HomogeneousDrag::Hf() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return
        tmp<volTensorField>
        (
            new volTensorField
            (
                IOobject
                (
                    "one",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedTensor("one", dimless, symmTensor::I)
            )
        );
}


// ************************************************************************* //

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


#include "dragClosure.H"
#include "phasePair.H"
#include "HDragCorrection.H"
#include "twoPhaseSystem.H"
#include "microDragClosure.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::dimensionSet Foam::dragClosure::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragClosure::dragClosure
(
    const bool registerObject,
    const phasePair& pair
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair)
{}


Foam::dragClosure::dragClosure
(
    const dictionary& dict,
    const phasePair& pair
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    pair_(pair),
    microDrag_
    (
        microDragClosure::New
        (
            dict.subDict("microscopicDragLaw"),
            pair
        )
    ),
    HCorrection_
    (
        HDragCorrection::New
        (
            dict.subDict("heterogeneousCorrection"),
            pair
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragClosure::~dragClosure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragClosure::writeFields() const
{
  volTensorField Hfact
  (
    IOobject
    (
      "drag.H",
      pair_.continuous().mesh().time().timeName(),
      pair_.continuous().mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    HCorrection_->Hf()
  );

  Hfact.write();

  volScalarField bMic
  (
    IOobject
    (
      "drag.betaMicro",
      pair_.continuous().mesh().time().timeName(),
      pair_.continuous().mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
              0.75                // it's 3/4
             *microDrag_->CdRe()
             *pair_.continuous().rho()
             *pair_.continuous().nu()
             /sqr(pair_.dispersed().d())
  );

  bMic.write();

}

Foam::tmp<Foam::volTensorField> Foam::dragClosure::Ki() const
{
    return
        0.75                // it's 3/4
       *microDrag_->CdRe()
       *HCorrection_->Hf()
       *pair_.continuous().rho()
       *pair_.continuous().nu()
       /sqr(pair_.dispersed().d());
}


Foam::tmp<Foam::volTensorField> Foam::dragClosure::K() const
{
    return max(pair_.dispersed(), pair_.dispersed().residualAlpha())*Ki();
}


Foam::tmp<Foam::surfaceTensorField> Foam::dragClosure::Kf() const
{
    return
        max
        (
            fvc::interpolate(pair_.dispersed()),
            pair_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki());
}

Foam::tmp<Foam::volScalarField> Foam::dragClosure::CdRe() const
{
  return microDrag_->CdRe();
}

bool Foam::dragClosure::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //

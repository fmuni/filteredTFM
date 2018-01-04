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

#include "stressClosure.H"

using namespace Foam;

//-------------------------- Constructors ---------------------------------//
Foam::StressClosure::StressClosure(  phaseModel&              phase,
                                     dictionary         closureDict,
                                     bool               isDispersed
                                 )
:
 phase_(phase),
 settings_(closureDict),
 mesoScaleStressClosure_
 (
   StressSubClosure::New
   (
    settings_.subDict("mesoScale"),
    phase_,
    word(settings_.subDict("mesoScale").lookup("type"))
   )
 ),
 microScaleStressClosure_
 (
   StressSubClosure::New
   (
    settings_.subDict("microScale"),
    phase_,
    word(settings_.subDict("microScale").lookup("type"))
   )
 ),
 isDispersed_(isDispersed)
{
  //Activate frictional just for the particle phase
  if(isDispersed_)
  {
    frictionalStressClosure_ =
      StressSubClosure::New
      (
       settings_.subDict("frictional"),
       phase_,
       word(settings_.subDict("frictional").lookup("type"))
      );
  }
};
//-------------------------- Destructors ----------------------------------//
Foam::StressClosure::~StressClosure()
{
};
//---------------------------   Methods  ----------------------------------//
void Foam::StressClosure::correct()
{
      mesoScaleStressClosure_->correct();
      microScaleStressClosure_->correct();

      if(frictionalStressClosure_.valid())
           frictionalStressClosure_->correct();

}

tmp<volScalarField> Foam::StressClosure::nuEff() const
{
    tmp<volScalarField> tmpNuEff
    (
        new volScalarField
        (
            IOobject
            (
                "nuEff."+phase_.name(),
                phase_.mesh().time().timeName(),
                phase_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase_.mesh(),
            dimensionedScalar("0",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
        )
    );

    volScalarField& nuEff = tmpNuEff.ref();

    mesoScaleStressClosure_->updateNu(nuEff);
    microScaleStressClosure_->updateNu(nuEff);

    if(isDispersed_)
     frictionalStressClosure_->updateNu(nuEff);

    return tmpNuEff;
}

tmp<volScalarField> Foam::StressClosure::k() const
{
  tmp<volScalarField> tmp (  mesoScaleStressClosure_->k()
                           + microScaleStressClosure_->k()
                          );
  if(isDispersed_)
   tmp = tmp + frictionalStressClosure_->k();

  return tmp;
}

tmp<volScalarField> Foam::StressClosure::muEff() const
{
 return  phase_.rho()*nuEff();
}


tmp<volScalarField> Foam::StressClosure::epsilon() const
{
  tmp<volScalarField> tmp (  mesoScaleStressClosure_->epsilon()
                           + microScaleStressClosure_->epsilon()
                          );
  if(isDispersed_)
   tmp = tmp + frictionalStressClosure_->epsilon();

  return tmp;
}

tmp<volSymmTensorField> Foam::StressClosure::R() const
{
  tmp<volSymmTensorField> tmp (  mesoScaleStressClosure_->R()
                               + microScaleStressClosure_->R()
                          );
  if(isDispersed_)
   tmp = tmp + frictionalStressClosure_->R();

  return tmp;
}

tmp<volScalarField> Foam::StressClosure::pPrime() const
{
    tmp<volScalarField> tmpPPrime
    (
        new volScalarField
        (
            IOobject
            (
                "nuEff."+phase_.name(),
                phase_.mesh().time().timeName(),
                phase_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase_.U().mesh().lookupObject<volScalarField>("p")*0
       )
    );

    volScalarField& pPrime = tmpPPrime.ref();

    mesoScaleStressClosure_->updateP(pPrime );
    microScaleStressClosure_->updateP(pPrime );

    if(isDispersed_)
     frictionalStressClosure_->updateP(pPrime );

    return tmpPPrime;
}

tmp<volSymmTensorField> Foam::StressClosure::devRhoReff() const
{
  tmp<volSymmTensorField> tmp (  mesoScaleStressClosure_->devRhoReff()
                               + microScaleStressClosure_->devRhoReff()
                          );
  if(isDispersed_)
   tmp = tmp + frictionalStressClosure_->devRhoReff();

  return tmp;
}

tmp<fvVectorMatrix> Foam::StressClosure::divDevRhoReff(volVectorField& U) const
{
  tmp<fvVectorMatrix> tmp (  mesoScaleStressClosure_->divDevRhoReff(U)
                               + microScaleStressClosure_->divDevRhoReff(U)
                          );
  if(isDispersed_)
   tmp = tmp + frictionalStressClosure_->divDevRhoReff(U);

  return tmp;
}

void Foam::StressClosure::writeFields() const
{
  mesoScaleStressClosure_->writeFields();
  microScaleStressClosure_->writeFields();

  if(isDispersed_)
   frictionalStressClosure_->writeFields();
}

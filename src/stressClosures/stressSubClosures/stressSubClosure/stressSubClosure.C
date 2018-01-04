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

#include "stressSubClosure.H"
#include "twoPhaseSystem.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "dimensionSet.H"

namespace Foam
{

 defineTypeNameAndDebug(StressSubClosure, 0);

 defineRunTimeSelectionTable(StressSubClosure, dictionary);
//-------------------------- Constructors ---------------------------------//
 StressSubClosure::StressSubClosure(
                                  const dictionary&          dict,
                                  phaseModel&               phase
                                 )
 :
 phase_(phase),
 settings_(dict),
 phaseThreshold_(settings_.lookupOrDefault("dispersedPhaseThreshold",0.005)),
 dispersedPhase_(phase.fluid().phase1())
 {
 };
//-------------------------- Destructors ----------------------------------//
 StressSubClosure::~StressSubClosure()
 {
 };
//---------------------------   Methods  ----------------------------------//

//Memory allocation methods
void StressSubClosure::createViscosity(word closureName)
{
  nuPrimePtr_.set
  (
    new volScalarField
    (
        IOobject
        (
            "nuPrime."+phase_.name() + "-" + closureName,
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("0",
                           dimensionSet(0,2,-1,0,0,0,0),
                           0.0
                          )
    )

  );
}

void StressSubClosure::createPressure(word closureName)
{
  pPrimePtr_.set
  (
    new volScalarField
    (
        IOobject
        (
            "pPrime."+phase_.name()+ "-" + closureName,
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_.U().mesh().lookupObject<volScalarField>("p")*0
    )

  );
}

void StressSubClosure::createLambda(word closureName)
{

     lambdaPrimePtr_.set
     (
       new volScalarField
       (
           IOobject
           (
               "lambdaPrime."+phase_.name()+ "-" + closureName,
               phase_.U().mesh().time().timeName(),
               phase_.U().mesh(),
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           phase_.U().mesh(),
           dimensionedScalar("0",
                              dimensionSet(0,2,-1,0,0,0,0),
                              0.0
                             )
       )

     );

}

void StressSubClosure::createASigma(word closureName)
{
  aSigmaPtr_.set
  (
    new volVectorField
    (
        IOobject
        (
            "aSigma."+phase_.name()+ "-" + closureName,
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedVector("0",
                                 dimensionSet(0,2,-2,0,0,0,0),
                               vector::zero
                             )
    )

  );
}

void StressSubClosure::checkAutoPtr(word property, bool empty) const
{
  if(empty)
  {
    FatalError << property << " was not allocated in "
             << "StressSubClosure " << typeName_()
             << " for phase " << phase_.name() << "\n"
             << "It seems like the class was not coded properly...";
  }
}


const volScalarField& StressSubClosure::nuPrime() const
{
    checkAutoPtr("Viscosity",nuPrimePtr_.empty());
    return nuPrimePtr_();
}

volScalarField& StressSubClosure::nuPrime()
{
    checkAutoPtr("Viscosity",nuPrimePtr_.empty());
    return nuPrimePtr_();
}

const volScalarField&  StressSubClosure::lambdaPrime() const
{
    checkAutoPtr("Lambda",lambdaPrimePtr_.empty());
    return lambdaPrimePtr_();
}

volScalarField&  StressSubClosure::lambdaPrime()
{
    checkAutoPtr("Lambda",lambdaPrimePtr_.empty());
    return lambdaPrimePtr_();
}

const volVectorField&  StressSubClosure::aSigma() const
{
    checkAutoPtr("Anisotropic Component",aSigmaPtr_.empty());
    return aSigmaPtr_();
}

volVectorField&  StressSubClosure::aSigma()
{
     checkAutoPtr("Anisotropic Component",aSigmaPtr_.empty());
    return aSigmaPtr_();
}

const volScalarField& StressSubClosure::pPrime() const
{
    checkAutoPtr("Pressure",pPrimePtr_.empty());
    return pPrimePtr_();
}

volScalarField& StressSubClosure::pPrime()
{
    checkAutoPtr("Pressure",pPrimePtr_.empty());
    return pPrimePtr_();
}

Foam::tmp<Foam::volScalarField>
StressSubClosure::k() const
{
    NotImplemented;
    return nuPrime();
}


Foam::tmp<Foam::volScalarField>
StressSubClosure::epsilon() const
{
    NotImplemented;
    return nuPrime();
}

Foam::tmp<Foam::volSymmTensorField>
StressSubClosure::R() const
{
  return R(phase_.U());
}

Foam::tmp<Foam::volSymmTensorField>
StressSubClosure::R(volVectorField& U) const
{
  const surfaceScalarField&   phi = phase_.phi();

    tmp<volSymmTensorField> tmpRt
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R."+phase_.name()+ "-" + typeName_(),
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase_.U().mesh(),
            dimensionedSymmTensor("R",
                                   dimensionSet(0,2,-2,0,0,0,0),
                                   symmTensor::zero
                                 )
        )
    );

    volSymmTensorField& Rt = tmpRt.ref();

    //Add viscous component to tensor if valid viscosity
    if(nuPrimePtr_.valid())
    {
     Rt -=  (nuPrime())*dev(twoSymm(fvc::grad(U)));
  ;  }

    //Add compressible component to tensor if valid compressible viscosity
    if(lambdaPrimePtr_.valid())
    {
     Rt -=  (lambdaPrime()*fvc::div(phi))*symmTensor::I;
    }

    //Add anisotropic component if available
    if(aSigmaPtr_.valid())
    {
     Rt +=  symm(aSigma()*vector::one);
    }

    return tmpRt;
}


Foam::tmp<Foam::surfaceScalarField>
StressSubClosure::pPrimef() const
{
    return fvc::interpolate(pPrime());
}

Foam::tmp<Foam::volSymmTensorField>
StressSubClosure::devRhoReff() const
{
    return devRhoReff(phase_.U());
}

Foam::tmp<Foam::volSymmTensorField>
StressSubClosure::devRhoReff(volVectorField& U) const
{

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", U.group()),
                phase_.mesh().time().timeName(),
                phase_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase_.rho()*R(U)()
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
StressSubClosure::divDevRhoReff
(
    volVectorField& U
) const
{

    Foam::tmp<Foam::fvVectorMatrix> tmpEq
    (
      new fvVectorMatrix(U,dimensionSet(1,1,-2,0,0,0,0))
    );

    fvVectorMatrix& Eq = tmpEq.ref();

    Eq += fvc::div(devRhoReff(U)());

    //Add laplacian if valid viscosity
    if(nuPrimePtr_.valid())
    {
     Eq -= fvm::laplacian(phase_.rho()*nuPrime(), U);
    }

    return tmpEq;
}

void StressSubClosure::updateNu(volScalarField& nu) const
{
    if(nuPrimePtr_.valid())
     nu += nuPrime();
}

void StressSubClosure::updateP(volScalarField& p) const
{
    if(pPrimePtr_.valid())
     p += pPrime();
}

void StressSubClosure::updatePf(surfaceScalarField& pf) const
{
    if(pPrimePtr_.valid())
     pf += pPrimef();
}

void StressSubClosure::updateLambda(volScalarField& lambda) const
{
    if(lambdaPrimePtr_.valid())
     lambda += lambdaPrime();
}

void StressSubClosure::updateASigma(volVectorField& ASigma) const
{
    if(aSigmaPtr_.valid())
     ASigma += aSigma();
}

const DynamicParameters& StressSubClosure::markers() const
{
  return phase_.fluid().dynPar();
}

void StressSubClosure::writeFields() const
{
  if(nuPrimePtr_.valid())
   nuPrime().write();

  if(lambdaPrimePtr_.valid())
   lambdaPrime().write();

  if(pPrimePtr_.valid())
   pPrime().write();

  if(aSigmaPtr_.valid())
   aSigma().write();
}


}

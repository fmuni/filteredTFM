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

#include "BoyerFrictional.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "error.H"

#define TOL 1e-10

namespace Foam
{

 defineTypeNameAndDebug(BoyerFrictional, 0);

 addToRunTimeSelectionTable(StressSubClosure, BoyerFrictional, dictionary);
//-------------------------- Constructors ---------------------------------//
 BoyerFrictional::BoyerFrictional(
                                  const dictionary&          dict,
                                  phaseModel&               phase
                                 )
 :
 StressSubClosure(dict,phase),
 mu1_(0.32),
 mu2_(0.7),
 I0_(5e-03),
 phiMax_(0.585),
 aLambda_(1.05,1,0.65)
 {
    if(phase_.isDispersed())
    {
        createViscosity(typeName_());
        createASigma(typeName_());
    }
    else
    {
      FatalError << "BoyerFrictional stress model should be used "
                 << "for the dispersed phase, not the continuous one!"
                 << exit(FatalError);
    }
 };
//-------------------------- Destructors ----------------------------------//
 BoyerFrictional::~BoyerFrictional()
 {

 };
//---------------------------   Methods  ----------------------------------//
 void BoyerFrictional::correct()
 {

     volScalarField        A(phiMax_/(phase_ + TOL) - 1.0 );
     volScalarField        invA(1.0/A);
     volScalarField        muf(phase_.fluid().otherPhase(phase_).mu());
     volSymmTensorField    relRefAT(markers().relAniTensor(phase_,aLambda_));


     nu() = muf / phase_.rho() *
              (
                   invA* 5.0/2.0 * phiMax_
                +  invA*invA*(mu1_ + (mu2_ - mu1_)/(1.0+I0_*invA*invA ))
              );

     aSigma() = relRefAT *
             (
               muf / phase_.rho() * markers().strainRate(phase_) * invA * invA
             ) ;

 }
}

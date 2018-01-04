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

#include "CloeteFrictional.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "error.H"

#define TOL 1e-10

namespace Foam
{

 defineTypeNameAndDebug(CloeteFrictional, 0);

 addToRunTimeSelectionTable(StressSubClosure, CloeteFrictional, dictionary);
//-------------------------- Constructors ---------------------------------//
 CloeteFrictional::CloeteFrictional(
                                  const dictionary&          dict,
                                  phaseModel&               phase
                                 )
 :
 StressSubClosure(dict,phase),
 x1_(4.016),
 x2_(-0.005543),
 x3_(0.1905),
 x4_(1.939),
 x5_(1.658),
 x6_(0.03935),
 x7_(12.78),
 x8_(5.084),
 DFref_(0.1286)
 {
  //Only frictional pressure
  createPressure(typeName_());
 };
//-------------------------- Destructors ----------------------------------//
 CloeteFrictional::~CloeteFrictional()
 {

 };
//---------------------------   Methods  ----------------------------------//
 void CloeteFrictional::correct()
 {

  scalarField          SrStar(
                               markers().strainRate(phase_)
                             * markers().settlingU()
                             / markers().gValue()
                            );

  scalar scaleP( markers().settlingU()*markers().settlingU() );

  forAll(SrStar,celli)
  {
    pPrime()[celli] =   scaleP*phase_.rho()[celli]
                      * pow(phase_[celli] + TOL,x1_)
                      * pow( (markers().filterSize(celli) - DFref_), x2_)
                      *(
                           x3_*pow(SrStar[celli]+ TOL,x4_)
                         * pow( (markers().filterSize(celli) - DFref_), x5_)
                         + x6_* pow(phase_[celli]+ TOL, x7_)
                         / pow( phase_.alphaMax() - phase_[celli]+ TOL, x8_)
                      );

  }

 }
}

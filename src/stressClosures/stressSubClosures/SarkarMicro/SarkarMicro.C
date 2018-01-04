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

#include "SarkarMicro.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"

#define TOL 1E-10



namespace Foam
{

 defineTypeNameAndDebug(SarkarMicro, 0);

 addToRunTimeSelectionTable(StressSubClosure, SarkarMicro, dictionary);
//-------------------------- Constructors ---------------------------------//
 SarkarMicro::SarkarMicro(
                                  const dictionary&          dict,
                                  phaseModel&               phase
                                 )
 :
 StressSubClosure(dict,phase)
 {
   createViscosity(typeName_());
   createPressure(typeName_());

 };
//-------------------------- Destructors ----------------------------------//
 SarkarMicro::~SarkarMicro()
 {
 };
//---------------------------   Methods  ----------------------------------//
 void SarkarMicro::correct()
 {

    //get strain rate Sr
    scalarField Sr(  markers().strainRate(phase_) );

    //Get gravity
    const scalar g(mag(phase_.fluid().g()).value());

    //Get settling velocity
    const scalar settlingU(phase_.fluid().dynPar().settlingU());

    scalar Cnu(
               0.00307
              *pow(
                    g/(settlingU*settlingU),
                   -6.0/7.0
                  )
              );

    forAll(Sr,celli)
    {
     //bounded delta alpha (ensures stability)
      scalar deltaAlpha(
                          max(
                             phase_.alphaMax() - phase_[celli],
                           phase_.residualAlpha().value()
                        )
     );


      //evaluate pressure
        pPrime()[celli] =  phase_.rho()()[celli]*0.01797
                        * markers().filterSize(celli,DIMENSIONAL)
                        * markers().filterSize(celli,DIMENSIONAL)
                        * Sr[celli]*Sr[celli]
                        * pow(max(phase_[celli],0.),1.645)
                        / deltaAlpha;

       //evaluate viscosity
       nuPrime()[celli] =  Cnu
                         * pow(markers().filterSize(celli,DIMENSIONAL),8.0/7.0)
                         * Sr[celli]
                         * pow(max(phase_[celli],0.),1.544)
                         / deltaAlpha;

    }

       //Apply smoothing to zero in case of dilute particle phase
       smoothToZero(pPrime());
       smoothToZero(nuPrime());


 }
}

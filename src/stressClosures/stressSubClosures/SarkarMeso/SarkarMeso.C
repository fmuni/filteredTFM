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

#include "SarkarMeso.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "error.H"
#include "stressSubClosure.H"


#define TOL 1E-10

namespace Foam
{

 defineTypeNameAndDebug(SarkarMeso, 0);

 addToRunTimeSelectionTable(StressSubClosure, SarkarMeso, dictionary);
//-------------------------- Constructors ---------------------------------//
 SarkarMeso::SarkarMeso(
                                  const dictionary&          dict,
                                  phaseModel&               phase
                                 )
 :
 StressSubClosure(dict,phase)
 {
   if(phase.isDispersed())
   {
    createViscosity(typeName_());
    createPressure(typeName_());
    createASigma(typeName_());
   }
   else
   {
    createViscosity(typeName_());
    createPressure(typeName_());
   }

 };
//-------------------------- Destructors ----------------------------------//
 SarkarMeso::~SarkarMeso()
 {
 };
//---------------------------   Methods  ----------------------------------//
 void SarkarMeso::correct()
 {
   if(phase_.isDispersed()) correctDispersed();
   else correctContinuous();
 }

 void SarkarMeso::correctDispersed()
 {

    //get strain rate Sr
    scalarField Sr(  markers().strainRate(phase_) );
    //Get gravity
    const scalar g(markers().gValue());
    //Get settling velocity
    const scalar settlingU(markers().settlingU());

    scalar Cnu(
               0.02518
              *pow(
                    g/(settlingU*settlingU),
                   -2./7.
                  )
              );


    scalar Cp(
               0.0236
              *pow(
                    g/(settlingU*settlingU),
                   3./7.
                  )
              );


    forAll(Sr,celli)
    {

 /*     //Set to zero if phase is less than tolerance
      if(phase_[celli] <  phase_.residualAlpha().value() )
      {
         pPrime()[celli] = pPrime()[celli]*0.;
         nuPrime()[celli] =  nuPrime()[celli]*0.;
         aSigma()[celli] = aSigma()[celli]*0.;
         continue;

      }
*/
      //bounded delta alpha (ensures stability)
      scalar deltaAlpha(
                          max(
                              phase_.alphaMax() - phase_[celli],
                              phase_.residualAlpha().value()
                           )
      );


     //evaluate pressure
      pPrime()[celli] =  phase_.rho()[celli]*Cp
                        * pow(markers().filterSize(celli,DIMENSIONAL),17./7.)
                        * Sr[celli]*Sr[celli]
                        * pow(fabs(phase_[celli]) +TOL,1.115)
                        / deltaAlpha;

      //evaluate viscosity
      nuPrime()[celli] =  Cnu
                        * pow(markers().filterSize(celli,DIMENSIONAL),12./7.)
                        * Sr[celli]
                        * pow(fabs(phase_[celli])+TOL,1.123)
                        / deltaAlpha;

     //evaluate anisotropic component
      aSigma()[celli].xx() =  -1./3. * 2.6 * pPrime()[celli]
                              / phase_.rho()[celli]
                              * pow (
                                      1.0 - min(phase_[celli] / phase_.alphaMax(),1.0) +TOL ,
                                      1.2
                                    );

      aSigma()[celli].zz() = aSigma()[celli].xx();
      aSigma()[celli].yy() = -2.*aSigma()[celli].xx();

    }

   //Apply smoothing to zero in case of dilute particle phase
   smoothToZero(pPrime());
   smoothToZero(nuPrime());
   smoothToZero(aSigma());


 }

 void SarkarMeso::correctContinuous()
 {


     //get strain rate Sr
     const scalarField& Sr(  markers().strainRate(phase_) );
    //Get gravity
    const scalar g(markers().gValue());

    //Get settling velocity
    const scalar settlingU(markers().settlingU());

    scalar Cp(
               pow(
                    g/(settlingU*settlingU),
                   5./7.
                  )
              );


    forAll(Sr,celli)
    {
      //Set to zero if only fluid
 /*    if(phase_[celli] > ( scalar(1) - phase_.fluid().otherPhase(phase_).residualAlpha().value() ) )
     {
       pPrime()[celli] = pPrime()[celli]*0.;
       nuPrime()[celli] =  nuPrime()[celli]*0.;
       continue;

     }
*/

      scalar phiSolid(1.-phase_[celli]);
      //evaluate pressure
      pPrime()[celli] =  phase_.rho()()[celli]*Cp
                        * pow(markers().filterSize(celli,DIMENSIONAL),19./7.)
                        * Sr[celli]*Sr[celli]
                        * (
                              0.0661
                            + 0.0164 * phiSolid
                            - 0.194 * phiSolid * phiSolid
                          );

      //evaluate viscosity
      nuPrime()[celli] =  pow(markers().filterSize(celli,DIMENSIONAL),2)
                        * Sr[celli]
                        * (
                              0.0330
                            + 0.218 * phiSolid
                            - 0.485 * phiSolid * phiSolid
                          );


    }

    smoothToZero(pPrime());
    smoothToZero(nuPrime());

 }



}

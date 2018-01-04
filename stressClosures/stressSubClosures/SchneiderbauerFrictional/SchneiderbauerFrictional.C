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

#include "SchneiderbauerFrictional.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "error.H"

#define TOL 1E-12
#define MAXVAULE 1e12

namespace Foam
{

 defineTypeNameAndDebug(SchneiderbauerFrictional, 0);

 addToRunTimeSelectionTable(StressSubClosure, SchneiderbauerFrictional, dictionary);
//-------------------------- Constructors ---------------------------------//
 SchneiderbauerFrictional::SchneiderbauerFrictional(
                                  const dictionary&          dict,
                                  phaseModel&               phase
                                 )
 :
 StressSubClosure(dict,phase),
 I0_(dict.lookupOrDefault("I0",0.279)), // see text below Eqn. 24 in Schneiderbauer et al.
 muSt_(dict.lookupOrDefault("muSt",0.38186)), //tan(20.9 degree), see text below Eqn. 24 in Schneiderbauer et al.
 muC_(dict.lookupOrDefault("muC",  0.64347)), //tan(32.76 degree), see text below Eqn. 24 in Schneiderbauer et al.
 muMax_(dict.lookupOrDefault("muMax",1e+04)), //see Table 3 of Schneiderbauer et al.
 pMax_(dict.lookupOrDefault("pMax",1e+04)), //see Table 3 of Schneiderbauer et al.
 b_(dict.lookupOrDefault("bJop",0.2)) //see text below Eqn. 31 in Schneiderbauer et al.
 {

   createViscosity(typeName_());
   createPressure(typeName_());

 };
//-------------------------- Destructors ----------------------------------//
 SchneiderbauerFrictional::~SchneiderbauerFrictional()
 {

 };
//---------------------------   Methods  ----------------------------------//
 void SchneiderbauerFrictional::correct()
 {
  scalarField   Sr( markers().strainRate(phase_) );

  scalarField  diam(phase_.d()());


  forAll(Sr,celli)
  {
    //bound shear rate to avoid incorrect behavior for extremely large shear rates
    Sr[celli] = min(Sr[celli],MAXVAULE);

    //regulartization of P_s_f calculation
    //Warning, there is a TYPO in Eqn. 33 of Schneiderbauer et al.!!
    scalar   X_reg = phase_.alphaMax()  - phase_[celli];
             X_reg = X_reg * X_reg + TOL; //ensure value is small and positive

    scalar   Y_reg = b_  * diam[celli] * Sr[celli];
             Y_reg = 4.0 * phase_.rho()[celli] * Y_reg * Y_reg + TOL; //ensure value is small and positive

    pPrime()[celli] = regularizeExp( X_reg, Y_reg, pMax_)
                    * Y_reg / X_reg;

    //regulartization of MU_s_f calculation
    //We use here Eqn. 32 inserted in Eqn. 27 of Schneiderbauer et al.
    //This means we evaluate: mu_s_fr = 2 mui rho_s  * (b * ds)^2 * sr / (phase_.alphaMax()  - phase_[celli])^2
    scalar Is(
                   2. * diam[celli] * Sr[celli] 
                  /
                  sqrt( pPrime()[celli] / phase_.rho()[celli] + TOL) //Need to add TOL here since pPrime can be zero!
    );

    //Avoid division by zery for small I_s! No regulartization needed here, since limit for small I_s is Mu_ist
    scalar mui(
                   muSt_
                   +
                   (muC_ - muSt_)
                   /
                   (
                      I0_/(Is+TOL) + 1.
                   )
    );
    //NOTE: X_reg can be reused from above!
    Y_reg =  b_ * diam[celli];
    Y_reg =  2.0 * mui * phase_.rho()[celli] * Y_reg * Y_reg * Sr[celli] + TOL;

    nuPrime()[celli] = regularizeExp(X_reg, Y_reg, muMax_)
                     * Y_reg / X_reg / phase_.rho()[celli]; //must return the KINEMATIC viscosity!
//    Info << "pPrime()[" << celli << "]: " << pPrime()[celli] 
//         << ", nuPrime()[" << celli << "]: " << nuPrime()[celli] 
//         << ", X_reg: " << X_reg
//         << ", Sr[celli]: " << Sr[celli]
//         << endl;
  }
}

scalar SchneiderbauerFrictional::regularizeExp(const scalar& argX, const scalar& argY, const scalar& argZ)
{
// This regularization is based on Schneiderbauer et al., Eq. 30
// Assumes that argY is NONZERO and POSITIVE
        return(
                scalar(1.0) - exp( -argZ*argX/argY )
              );
}
}

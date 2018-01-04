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

#include "KEquilibrium.H"
#include "twoPhaseSystem.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "error.H"
#include "dynamicParameters.H"
#include "addToRunTimeSelectionTable.H"
#include "auxEquations.H"

#define TOL 1e-20

using namespace Foam;


    defineTypeNameAndDebug(KEquilibrium, 0);
    addToRunTimeSelectionTable(EqnSubClosure, KEquilibrium, dictionary);
//-------------------------- Constructors ---------------------------------//
 KEquilibrium::KEquilibrium(
                                  twoPhaseSystem&      fluid,
                                  dictionary&          dict
                                 )
 :
 EqnSubClosure(fluid,dict),
 otherKName_(dict.lookup("otherPhaseK")),
 Cnug_(0.4),
 Cepg_(1.),
 Cnus_(0.25),
 Ceps_(1.),
 epgs_(0.8)
 {
   //It requires an equation for the kinetic energy of the other phase
   fluid_.equations().requestEquation(otherKName_,typeName_());
   word currPhase = word(dict.lookup("phase"));

   if(currPhase=="continuous")
   {
     phaseId_ = CONTINUOUS;
   }
   else if (currPhase=="dispersed")
   {
     phaseId_ = DISPERSED;
   }
   else
   {
     FatalError << "Invalid phase selected for KEquilibrium closure. \n"
                << "Valid phase entries for \"phase\" are \"dispersed\" or \"continuous\".";

   }
 };
//-------------------------- Destructors ----------------------------------//
 KEquilibrium::~KEquilibrium()
 {
 };
//---------------------------   Methods  ----------------------------------//
void KEquilibrium::closeEquation(fvMatrix<scalar>& Eqn, volScalarField& field)
{

 if(phaseId_==DISPERSED)
 {
   phaseModel& alpha =  fluid_.phase1();
   scalar C1 = Cnus_;
   scalar C2 = Ceps_;
   kClosure(Eqn,field,alpha,C1,C2);
 }
 else if(phaseId_==CONTINUOUS)
 {
   phaseModel& alpha =  fluid_.phase2();
   scalar C1 = Cnug_;
   scalar C2 = Cepg_;
   kClosure(Eqn,field,alpha,C1,C2);
 }

}

void KEquilibrium::kClosure(
                             fvMatrix<scalar>&     Eqn,
                             volScalarField&     field,
                             phaseModel&        alphaM,
                             scalar                 C1,
                             scalar                 C2
                           )
{
  //Get strainRate
  const volScalarField& Sr(markers().strainRate(alphaM));

  //Get average drag
  volScalarField  beta(fluid_.KdAve()());

  //Auxiliary relation
  scalar C2Square(C2*C2);

  //Get other phase K
  const volScalarField* kd(fluid_.equations().getFieldScalar(otherKName_));

  //Correct equation
  forAll(Sr,celli)
  {
    Eqn.diag()[celli] += C2Square;

    scalar lm(C1*markers().filterSize(celli));
    scalar alpha(min(alphaM[celli],TOL));
    scalar betaFac(- beta[celli]*lm/(alpha*alphaM.rho()[celli]));

    scalar source(    -betaFac
                      + sqrt(
                           betaFac*betaFac
                         + (lm*lm*Sr[celli]*Sr[celli])
                         + 2.*(
                             C2*betaFac*sqrt((*kd)[celli])

                         )

                       )
                 );
     Eqn.source()[celli] += source;
  }

}

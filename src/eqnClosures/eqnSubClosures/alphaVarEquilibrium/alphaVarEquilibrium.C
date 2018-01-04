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

#include "alphaVarEquilibrium.H"
#include "twoPhaseSystem.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "error.H"
#include "dynamicParameters.H"
#include "addToRunTimeSelectionTable.H"
#include "auxEquations.H"

#define TOL 1e-20

using namespace Foam;


    defineTypeNameAndDebug(alphaVarEquilibrium, 0);
    addToRunTimeSelectionTable(EqnSubClosure, alphaVarEquilibrium, dictionary);
//-------------------------- Constructors ---------------------------------//
 alphaVarEquilibrium::alphaVarEquilibrium(
                                  twoPhaseSystem&      fluid,
                                  dictionary&          dict
                                 )
 :
 EqnSubClosure(fluid,dict),
 KdName_(dict.lookup("dispersedPhaseK")),
 Cphis_(0.25),
 epphis_(0.1),
 Ceps_(1.0)
 {
   //It requires an equation for the kinetic energy of the continuous phase
   fluid_.equations().requestEquation(KdName_,typeName_());
 };
//-------------------------- Destructors ----------------------------------//
 alphaVarEquilibrium::~alphaVarEquilibrium()
 {
 };
//---------------------------   Methods  ----------------------------------//
void alphaVarEquilibrium::closeEquation(fvMatrix<scalar>& Eqn, volScalarField& field)
{
 scalarField divU(fvc::div(fluid_.phase1().phi()));
 const volScalarField* kd(fluid_.equations().getFieldScalar(KdName_));
 vectorField gradPhi(fvc::grad(fluid_.phase1()));


 forAll(fluid_.mesh().C(),celli)
 {
   Eqn.diag()[celli] += 3./8.;

   scalar lms(Cphis_*markers().filterSize(celli));

   scalar deriv(
               gradPhi[celli].x() + gradPhi[celli].y() + gradPhi[celli].z()
   );

   scalar kdFac(
           divU[celli]
         + Cphis_*Ceps_* sqrt((*kd)[celli])/lms
   );

   scalar source(    epphis_*epphis_*(*kd)[celli]
                     * deriv*deriv
                    /
                    (
                      kdFac*kdFac
                    )

                );
    Eqn.source()[celli] += source;
 }

}

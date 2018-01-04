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

#include "eqnSubClosure.H"
#include "twoPhaseSystem.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "error.H"
#include "dynamicParameters.H"

using namespace Foam;


 defineTypeNameAndDebug(EqnSubClosure, 0);

 defineRunTimeSelectionTable(EqnSubClosure, dictionary);
//-------------------------- Constructors ---------------------------------//
 EqnSubClosure::EqnSubClosure(
                                  twoPhaseSystem&            fluid,
                                  dictionary&          dict
                                 )
 :
 fluid_(fluid),
 settings_(dict)
 {
 };
//-------------------------- Destructors ----------------------------------//
 EqnSubClosure::~EqnSubClosure()
 {
 };
//---------------------------   Methods  ----------------------------------//
DynamicParameters&  EqnSubClosure::markers()
{
    return fluid_.dynPar();
}

void EqnSubClosure::closeEquation(fvMatrix<scalar>& Eqn, volScalarField& field)
{
  FatalError << "Closure " << typeName_()
             << " cannot close equation for scalar field " << field.name()
             << " because it was not designed to work on scalar equations!"
             << exit(FatalError);

}

void EqnSubClosure::closeEquation(fvMatrix<vector>& Eqn, volVectorField& field)
{
  FatalError << "Closure " << typeName_()
             << " cannot close equation for vector field " << field.name()
             << " because it was not designed to work on vector equations!"
             << exit(FatalError);

}

void EqnSubClosure::closeEquation(fvMatrix<tensor>& Eqn, volTensorField& field)
{
  FatalError << "Closure " << typeName_()
             << " cannot close equation for tensor field " << field.name()
             << " because it was not designed to work on tensor equations!"
             << exit(FatalError);

}

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

#ifndef alphaVarEquilibrium_SUBCLOSURE_H
#define alphaVarEquilibrium_SUBCLOSURE_H

#include "fvCFD.H"
#include "fvMatrices.H"
#include "eqnSubClosure.H"
#include "volFields.H"
#include "geometricOneField.H"
#include "twoPhaseSystem.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

 class alphaVarEquilibrium : public EqnSubClosure
 {
  protected:

  word KdName_;

  //Molde constants
  scalar Cphis_,epphis_,Ceps_;


  public:

       //- Runtime type information
    TypeName("alphaVarEquilibrium");

    // Declare runtime constructor selection table

    // Constructors

        alphaVarEquilibrium
        (
              twoPhaseSystem& fluid,
              dictionary& dict
        );

    // Destructors

       virtual ~alphaVarEquilibrium();

    // Selector


    //Methods

    //Update equation
     void closeEquation(fvMatrix<scalar>& Eqn, volScalarField& field);

 };

}

#endif

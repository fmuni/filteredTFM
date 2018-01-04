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

#ifndef equation_SUBCLOSURE_H
#define equation_SUBCLOSURE_H

#include "fvCFD.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "geometricOneField.H"
#include "twoPhaseSystem.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

 class EqnSubClosure
 {
  protected:

  twoPhaseSystem&       fluid_;

  dictionary               settings_;

  //Shorthand for markers
  DynamicParameters& markers();


  public:

       //- Runtime type information
    TypeName("EqnSubClosure");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            EqnSubClosure,
            dictionary,
            (
              twoPhaseSystem& fluid,
              dictionary& dict
            ),
            (fluid,dict)
        );


    // Constructors

        EqnSubClosure
        (
              twoPhaseSystem& fluid,
              dictionary& dict
        );

    // Destructors

       virtual ~EqnSubClosure();

    // Selector

       static autoPtr<EqnSubClosure> New
       (
           twoPhaseSystem& fluid,
           dictionary& dict,
           word& subClosureType
       );

    //Methods

    //Equation closures (now are differentiated for appropiate type)
    virtual void closeEquation(fvMatrix<scalar>& Eqn, volScalarField& field);
    virtual void closeEquation(fvMatrix<vector>& Eqn, volVectorField& field);
    virtual void closeEquation(fvMatrix<tensor>& Eqn, volTensorField& field);

 };

}

#endif

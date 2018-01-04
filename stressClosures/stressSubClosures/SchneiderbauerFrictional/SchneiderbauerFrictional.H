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

Description
    This stress model is based on Schneiderbauer et al., 
    Chem Eng Sci 80(2012), 279–292
\*---------------------------------------------------------------------------*/

#ifndef SchneiderbauerFrictional_H
#define SchneiderbauerFrictional_H

#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "phaseModel.H"
#include "stressSubClosure.H"

namespace Foam
{

 class SchneiderbauerFrictional
 :
    public StressSubClosure
 {
  protected:

      //Model coefficients
      scalar I0_,muSt_,muC_,muMax_,pMax_,b_;

      scalar regularizeExp(const scalar& arg1, const scalar& arg2, const scalar& arg3);


  public:

       //- Runtime type information
    TypeName("SchneiderbauerFrictional");

     // Constructors

        SchneiderbauerFrictional
        (
            const dictionary& dict,
            phaseModel& phase
        );

    // Destructors

     virtual ~SchneiderbauerFrictional();


    //Correct stress tensor
    virtual void correct();

 };

}

#endif

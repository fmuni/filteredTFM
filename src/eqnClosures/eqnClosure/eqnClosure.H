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

#ifndef equation_CLOSURE_H
#define equation_CLOSURE_H

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include <vector>
#include "eqnSubClosure.H"

namespace Foam
{

 class EqnSubClosure;

 class EqnClosure
 {
  protected:

  twoPhaseSystem&              fluid_;

  dictionary                      settings_;

  word                           fieldName_;

  wordList                 subClosureNames_;

  std::vector<EqnSubClosure*>  subClosures_;

  public:


    //Members

    //Methods

     //Construct from phaseModel object and dictionary
     EqnClosure
     (
      twoPhaseSystem&         fluid,
      dictionary         closureDict,
      word               fieldname
     );

     //Destructor
     ~EqnClosure();

     //Correct sub closures
     void correct() {};

     //Close governing equation
     template<class T>
     void   closeEquation( fvMatrix<T>& Eqn,
                           GeometricField<T,fvPatchField,volMesh>&  field
                         )
     {
      //Call all sub closures
      for(unsigned int subC=0;subC<subClosures_.size();subC++)
      {
          subClosures_[subC]->closeEquation(Eqn,field);
      }
     }


 };
}

#endif

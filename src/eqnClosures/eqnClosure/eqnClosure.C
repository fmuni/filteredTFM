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

#include "eqnClosure.H"
#include "eqnSubClosure.H"
#include "dynamicParameters.H"

using namespace Foam;

//-------------------------- Constructors ---------------------------------//
EqnClosure::EqnClosure(  twoPhaseSystem&          fluid,
                         dictionary         closureDict,
                         word                 fieldname
                      )
:
 fluid_(fluid),
 settings_(closureDict),
 fieldName_(fieldname),
 subClosureNames_(closureDict.lookup("subClosures"))
{
  forAll(subClosureNames_,subC)
  {
    EqnSubClosure * tmp(

      EqnSubClosure::New(
         fluid_,
         settings_,
         subClosureNames_[subC]
     ).ptr()
    );

   subClosures_.push_back(tmp);

  }
};
//-------------------------- Destructors ----------------------------------//
EqnClosure::~EqnClosure()
{
    for(unsigned int subC=0;subC<subClosures_.size();subC++)
     delete subClosures_[subC];
};
//---------------------------   Methods  ----------------------------------//

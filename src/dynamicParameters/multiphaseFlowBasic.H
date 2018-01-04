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

#ifndef multiphaseFlowBasic_H
#define multiphaseFlowBasic_H

#include "fvCFD.H"

namespace Foam
{
 class multiphaseFlowBasic
 {
  protected:

   bool      verbose_;

   //-List of implemented drag laws
   wordList dragLaws_;

   //-Beetstra drag
   double F_Beetstra
   (
           double,
           double
             );
   //-Wen-Yu drag
   double
   F_WenYu
    (
           double,
           double
             );
   //Koch-Hill drag
   double
   F_KochHill
    (
           double,
           double
             );
   //Set settling velocity and Reynolds
   void
   uSettlingSphereSuspension( double dP,
                              double rhoP,
                              double g,
                              double etaFluid,
                              double rhoFluid,
                              double phiP,
                              word    dragLaw,
                              double &u,
                              double &Re
                              );



  public:

      // Constructor
   multiphaseFlowBasic();

   // Destructor
   ~multiphaseFlowBasic();

   //variables
   struct settlingParams
   {
       double  u;              //the settling velocity of an isolated particle
       double  Re;             //the Re of an isolated particle
       double  Tref;           //u/g, the reference time scale
       double  Lref;           //u^2/g, the reference length scale
       double  Lchar;          //u^2/g*FrP^-0.5, a characteristic length scale
       double  Lchar2;          //u^2/g*FrP^-0.5, a characteristic length scale
       double  FrP;             //u^2/g/dp, the particle Froude number
       double  FrPf;             //u^2/g/Df, the filter Froude number
   };

   settlingParams settling;


   //Set all settling parameters
   void
   setupSettling(             double dP,
                              double rhoP,
                              double g,
                              double etaFluid,
                              double rhoFluid,
                              word    dragLaw,
                              double filterSize
                              );




//TODO: Write settling parameters to file
  void writeSettling() {};

 };


}
#endif

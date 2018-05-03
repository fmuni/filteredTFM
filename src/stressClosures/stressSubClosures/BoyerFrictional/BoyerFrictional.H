/*---------------------------------------------------------------------------*\
    ███████╗████████╗███████╗███╗   ███╗
    ██╔════╝╚══██╔══╝██╔════╝████╗ ████║       Filtered
    █████╗     ██║   █████╗  ██╔████╔██║       Two
    ██╔══╝     ██║   ██╔══╝  ██║╚██╔╝██║       Fluid
    ██║        ██║   ██║     ██║ ╚═╝ ██║       Model
    ╚═╝        ╚═╝   ╚═╝     ╚═╝     ╚═╝

-------------------------------------------------------------------------------
License
    This file is part of FTFM.

    FTFM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    FTFM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Contributors
    Federico Municchi, Purdue University, 2018
    
Description
    Frictional stress model proposed by Boyer et al. (2011).
    Here the anisotropic pressure is decomposed into hydrostatic (pPrime)
    and deviatoric (aSigma) components.
    Notice that for the hydrostatic pressure, only the derivative with respect
    to the phase field is required.
\*---------------------------------------------------------------------------*/

#ifndef BoyerFrictional_H
#define BoyerFrictional_H

#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "phaseModel.H"
#include "stressSubClosure.H"

namespace Foam
{

 class BoyerFrictional
 :
    public StressSubClosure
 {
  protected:

      //Model coefficients
      scalar mu1_,mu2_,I0_;

      //Maximum particle volume fraction
      scalar phiMax_;

      //Anisotropy coefficients
      vector aLambda_;




  public:

       //- Runtime type information
    TypeName("BoyerFrictional");

     // Constructors

        BoyerFrictional
        (
            const dictionary& dict,
            phaseModel& phase
        );

    // Destructors

     virtual ~BoyerFrictional();


    //Correct stress tensor
    virtual void correct();

 };

}

#endif

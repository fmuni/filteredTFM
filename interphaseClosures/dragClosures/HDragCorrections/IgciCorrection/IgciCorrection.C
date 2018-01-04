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

#include "IgciCorrection.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "dynamicParameters.H"

#define TOL 1e-10
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HDragCorrections
{
    defineTypeNameAndDebug(IgciCorrection, 0);
    addToRunTimeSelectionTable(HDragCorrection, IgciCorrection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HDragCorrections::IgciCorrection::IgciCorrection
(
    const dictionary& dict,
    const phasePair& pair
)
:
    HDragCorrection(dict, pair),
    settings_(dict),
    HIgci
    (
      IOobject
      (
          "HIgci",
           pair_.dispersed().fluid().mesh().time().timeName(),
           pair_.dispersed().fluid().mesh()
      ),
      pair_.dispersed().fluid().mesh(),
      dimensionedScalar("zero", dimless, scalar(0))
    ),
    operationAlpha(settings_.lookupOrDefault<word>("phiOperation","none"))

{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HDragCorrections::IgciCorrection::~IgciCorrection()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::HDragCorrections::IgciCorrection::Hf() const
{

    const fvMesh&  mesh = pair_.dispersed().fluid().mesh();
    //Calculate heterogeneity factor
    forAll(mesh.C(),celli)
    {
        //Get filtered Freud number
        scalar freudFPow =
           pow(markers().filterFr(celli),-1.3);

        scalar fFactor =   freudFPow
                         / ( 1.5 + freudFPow );

       //solid volume fraction of current celli
       scalar alpha(fabs(pair_.dispersed()[celli])+TOL);

       if(alpha<0.0012)
       {
          HIgci[celli] =  1. - fFactor *
                          2.7*pow(alpha,0.234);
       }
       else if (  alpha < 0.014 )
       {
           HIgci[celli] = 1. - fFactor *
                          ( 0.963 -0.019*pow(alpha,-0.445));

       }
       else if ( alpha < 0.25 )
       {
           HIgci[celli] = 1.  - fFactor *
                          (
                            0.868*exp(-0.38*alpha)
                           -0.176*exp(-119.2*alpha)
                          );
       }
       else if ( alpha < 0.455 )
       {
           HIgci[celli] =  1.  - fFactor *
                          (
                            -0.0000459*exp(19.75*alpha)
                            +0.852*exp(-0.268*alpha)
                          );
       }
       else if (  alpha < 0.59 )
       {
           HIgci[celli] =  1.  - fFactor *
                          (
                            (alpha-0.59) *
                            ( - 1051.*alpha*alpha*alpha
                              + 2203.*alpha*alpha
                              - 1054.*alpha
                              + 162.
                            )
                          );
       }
       else
       {
          HIgci[celli] =  1.;
       }

    }
   //Return tensor field
   return HIgci*tensor::I;
}



// ************************************************************************* //

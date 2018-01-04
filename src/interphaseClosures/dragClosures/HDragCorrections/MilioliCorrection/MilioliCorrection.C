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

#include "MilioliCorrection.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "dynamicParameters.H"
#include "error.H"

#define TOL 1e-10
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HDragCorrections
{
    defineTypeNameAndDebug(MilioliCorrection, 0);
    addToRunTimeSelectionTable(HDragCorrection, MilioliCorrection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HDragCorrections::MilioliCorrection::MilioliCorrection
(
    const dictionary& dict,
    const phasePair& pair
)
:
    HDragCorrection(dict, pair),
    settings_(dict),
    HMilioli
    (
      IOobject
      (
          "HMilioli",
           alphaD().fluid().mesh().time().timeName(),
           alphaD().fluid().mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE
      ),
      alphaD().fluid().mesh(),
      dimensionedTensor("zero", dimless, tensor::zero)
    ),
    operationAlpha(settings_.lookupOrDefault<word>("phiOperation","none"))

{
   //Check if model is valid in this mesh
   const fvMesh& mesh = alphaD().fluid().mesh();

   forAll(mesh.C(),celli)
   {
     if(
        markers().filterSize(celli) > 25.
     )
     {
       FatalError << "Milioli's heterogeneity correction is only valid"
                 << "for filters larger than 25 particle diameters!\n"
                 << "A filter size of" << markers().filterSize(celli)
                 << " particle diameters was found!";
     }
   }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HDragCorrections::MilioliCorrection::~MilioliCorrection()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::HDragCorrections::MilioliCorrection::Hf() const
{

    const fvMesh&           mesh = alphaD().fluid().mesh();
    const volVectorField&   markerU(markers().scaledSlipU(
                                                          pair_.continuous(),
                                                          pair_.dispersed()
                                                         )
                                   );
    //Calculate heterogeneity factor
    forAll(mesh.C(),celli)
    {
      scalar alpha(fabs(alphaD()[celli]));

      scalar henv1(0.);
      scalar henv2(0.8428 +0.6393*alpha - 0.6743*alpha*alpha);
      vector h1(vector::zero);
      vector finf(vector::zero);
      vector fext(vector::zero);
      vector hlim(vector::zero);


      if(alpha<0.1)
      {
        henv1 = ( 0.5643*(1.0+alpha) * pow(alpha,0.15) )
               / ( 0.5766*pow(alpha,0.3) +0.1997 );
      }
      else if(alpha < 0.54)
      {
        henv1 = henv2;
      }
      else if(alpha < 0.65)
      {
        henv1 =  (0.4099*pow((alphaD().alphaMax() + alphaD().residualAlpha().value()-alpha),0.25))
                / (pow(alpha+TOL,-0.25)-0.9281);
      }
      else
      {
        henv2 = 0.;
      }

      for(int comp=0;comp<3;comp++)
      {
        scalar magU(mag(markerU[celli].component(comp)));
        scalar fSize(markers().filterSize(celli,DIMENSIONLESS));

        if(fSize<1.542)
        {
           h1.component(comp) =
            alpha*(
               1.076 + 0.12 * magU - 0.02*( magU + 0.01)
            )
            + (0.084 + 0.09 * magU - 0.01/(0.1*magU +0.01));

        }
        else if(fSize<3.084)
        {
            h1.component(comp) =
             alpha*(
                1.268 -0.2 * magU - 0.14*( magU + 0.01)
             )
             + (0.385 + 0.09 * magU - 0.05/(0.2*magU +0.01));

        }
        else if(fSize<6.168)
        {
            h1.component(comp) =
            alpha*(
                  (0.18*magU + 0.1)/(0.14*magU+0.01)
               )
               + (0.9454   - 0.09/(0.2*magU +0.01));

        }
        else if(fSize<10.28)
        {
             h1.component(comp) =
             alpha*(
                   (0.05*magU + 0.3)/(0.4*magU+0.06)
                )
                + (0.9466   - 0.05/(0.11*magU +0.01));

        }
        else if(fSize<14.392)
        {
             h1.component(comp) =
             alpha*(
                   (1.3*magU + 2.2)/(5.2*magU+0.07)
                )
                + (0.9363   - 0.11/(0.3*magU +0.01));

        }
        else if(fSize<18.504)
        {
             h1.component(comp) =
             alpha*(
                   (2.6*magU + 4.)/(10.*magU+0.08)
                )
                + (0.926   - 0.17/(0.5*magU +0.01));

        }
        else
        {
             h1.component(comp) =
             alpha*(
                   (2.5*magU + 4.)/(10.*magU+0.08)
                )
                + (0.9261 +  - 0.17/(0.5*magU +0.01));

        }


         if(h1.component(comp)>0.)
        {
          finf.component(comp) =  0.882
                   * ( 2.145
                      - (
                          7.8* pow(magU  ,1.8)
                         /(0.5586 + 7.746*pow(magU ,1.8))
                        )
                     );
          fext.component(comp) =  finf.component(comp)
                                 *
                                 (
                                   3.494/(1. +8.*pow(fSize,0.4))
                                   +
                                   0.882
                                 );

          hlim.component(comp) = fext.component(comp)*h1.component(comp);
        }

      }


      HMilioli[celli].xx() = 1.0 - min (henv1,hlim.component(0));
      HMilioli[celli].yy() = 1.0 - min (henv1,hlim.component(1));
      HMilioli[celli].zz() = 1.0 - min (henv1,hlim.component(2));

    }
   //Return tensor field
   return HMilioli;
}



// ************************************************************************* //

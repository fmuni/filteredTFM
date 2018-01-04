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

#ifndef IgciCorrection_H
#define IgciCorrection_H

#include "HDragCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace HDragCorrections
{

/*---------------------------------------------------------------------------*\
                           Class IgciCorrection Declaration
\*---------------------------------------------------------------------------*/

class IgciCorrection
:
    public HDragCorrection
{
  protected:

  dictionary  settings_;
  
  //Igci's correction is isotropic
   mutable volScalarField    HIgci;

  //-Name oft the operation for alpha field
  // default=="none"
  word operationAlpha;


  public:

    //- Runtime type information
    TypeName("Igci");


    // Constructors

        //- Construct from a dictionary and a phase pair
        IgciCorrection
        (
            const dictionary& dict,
            const phasePair& pair
        );



    //- Destructor
    virtual ~IgciCorrection();


    // Member Functions

        //- Heterogeneous correction
        tmp<volTensorField> Hf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HDragCorrections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

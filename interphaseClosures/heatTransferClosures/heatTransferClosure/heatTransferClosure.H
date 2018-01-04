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

#ifndef heatTransferClosure_H
#define heatTransferClosure_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;
class microHeatTransferClosure;
class HHeatTransferCorrection;


/*---------------------------------------------------------------------------*\
                           Class heatTransferClosure Declaration
\*---------------------------------------------------------------------------*/

class heatTransferClosure
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //- microscopic drag law
        autoPtr<microHeatTransferClosure>    microHeat_;

        //- Heterogeneous drag correction
        autoPtr<HHeatTransferCorrection>           HCorrection_;


public:


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        // Construct without residual constants
        heatTransferClosure
        (
            const phasePair& pair,
            const bool registerObject
        );

        // Construct with residual constants
        heatTransferClosure
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    ~heatTransferClosure();




    // Member Functions


        //Return Heat transfer Coefficient
        tmp<volScalarField> K() const;


        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

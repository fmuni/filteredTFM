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

#ifndef MICRO_DRAG_CLOSURE_H
#define MICRO_DRAG_CLOSURE_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;


/*---------------------------------------------------------------------------*\
                           Class microDragClosure Declaration
\*---------------------------------------------------------------------------*/

class microDragClosure
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("microDragClosure");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            microDragClosure,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        microDragClosure
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~microDragClosure();


    // Selectors

        static autoPtr<microDragClosure> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const = 0;


        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

#ifndef HDragCorrection_H
#define HDragCorrection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
#include "dynamicParameters.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                           Class HDragCorrection Declaration
\*---------------------------------------------------------------------------*/

class HDragCorrection
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //Shorthand for dispersed phase fraction
        const phaseModel&  alphaD() const;

        //Shorthand for continuous phase fraction
        const phaseModel&  alphaC() const;

        //Shorthand for dynamic parameters
        const DynamicParameters& markers() const;


public:

    //- Runtime type information
    TypeName("HDragCorrection");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        HDragCorrection,
        dictionary,
        (
            const dictionary& dict,
            const phasePair& pair
        ),
        (dict, pair)
    );


    // Constructors

        //- Construct from a dictionary and a phase pair
        HDragCorrection
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~HDragCorrection();


    // Selectors

        static autoPtr<HDragCorrection> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Heterogeneous factor
        virtual tmp<volTensorField> Hf() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

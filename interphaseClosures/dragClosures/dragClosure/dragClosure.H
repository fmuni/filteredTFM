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

#ifndef DRAG_CLOSURE_H
#define DRAG_CLOSURE_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"


namespace Foam
{

class phasePair;
class microDragClosure;
class HDragCorrection;


/*---------------------------------------------------------------------------*\
                           Class dragClosure Declaration
\*---------------------------------------------------------------------------*/

class dragClosure
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //- microscopic drag law
        autoPtr<microDragClosure>    microDrag_;

        //- Heterogeneous drag correction
        autoPtr<HDragCorrection>   HCorrection_;


public:


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        // Construct without residual constants
        dragClosure
        (
            const bool registerObject,
            const phasePair& pair
        );

        // Construct with residual constants
        dragClosure
        (
            const dictionary& dict,
            const phasePair& pair
        );

    //- Destructor
    ~dragClosure();




    // Member Functions

        //- Write drag tensor
        void writeFields() const;

        //- Return the phase-intensive drag coefficient Ki
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(U2-U1)
        tmp<volTensorField> Ki() const;

        //- Return the drag coefficient K
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(U2-U1)
        tmp<volTensorField> K() const;

        //- Return the drag coefficient Kf
        //  used in the face-momentum equations
         tmp<surfaceTensorField> Kf() const;

        //- Return microscopic drag coefficient
        tmp<volScalarField>    CdRe() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

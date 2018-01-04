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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::twoPhaseSystem::mesh() const
{
    return mesh_;
}

inline const Foam::dimensionedVector& Foam::twoPhaseSystem::g() const
{
    return g_;
}


inline const Foam::phaseModel& Foam::twoPhaseSystem::phase1() const
{
    return phase1_;
}


inline Foam::phaseModel& Foam::twoPhaseSystem::phase1()
{
    return phase1_;
}


inline const Foam::phaseModel& Foam::twoPhaseSystem::phase2() const
{
    return phase2_;
}


inline Foam::phaseModel& Foam::twoPhaseSystem::phase2()
{
    return phase2_;
}


inline const Foam::phaseModel& Foam::twoPhaseSystem::otherPhase
(
    const phaseModel& phase
) const
{
    if (&phase == &phase1_)
    {
        return phase2_;
    }
    else
    {
        return phase1_;
    }
}


inline const Foam::surfaceScalarField& Foam::twoPhaseSystem::phi() const
{
    return phi_;
}


inline Foam::surfaceScalarField& Foam::twoPhaseSystem::phi()
{
    return phi_;
}


inline const Foam::volScalarField& Foam::twoPhaseSystem::dgdt() const
{
    return dgdt_;
}


inline Foam::volScalarField& Foam::twoPhaseSystem::dgdt()
{
    return dgdt_;
}


inline Foam::tmp<Foam::surfaceScalarField>& Foam::twoPhaseSystem::pPrimeByA()
{
    return pPrimeByA_;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling
    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
                                Copyright (C) 2017-     Graz University of  
                                                        Technology, IPPT
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.
    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    eulerianFilteredTFM

Description
    Solver for a system of 2 compressible fluid phases with one phase
    dispersed, e.g., solid particles in a gas. The solver invovles latest
    closures for the so-called "filtered TFM" approach, which accounts for 
    the effect of unresolved mesoscale structures. Thus, the present solver 
    can be used on comparably coarse computational grids without discriminating
    significantly the fidelity of the simulation result.
    
    The solver is based on, but not identical to, OpenFOAM's twoFluidEulerFOAM solver.
    
Contributors
    Federico Municchi, TU Graz, 2017
    Stefan Radl, TU Graz, 2017 -

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    Switch implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).lookupOrDefault<Switch>
        (
            "implicitPhasePressure", false
        )
    );

    #include "pUf/createDDtU.H"
    #include "pU/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            #include "contErrs.H"

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
              //  #include "EEqns.H"
                #include "pUf/pEqn.H"
                #include "pUf/DDtU.H"
            }
            else
            {
                #include "pU/UEqns.H"
              //  #include "EEqns.H"
                #include "pU/pEqn.H"
                #include "pU/DDtU.H"
            }

   //         if (pimple.turbCorr())
   //         {
                fluid.correctStressClosures();
     //       }
        }

        #include "write.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey SIMURZIN
##


#---------------------------------------------------------------------------
from Foam import ref, man


#---------------------------------------------------------------------------
def _createFields( runTime, mesh ):
    
    ref.ext_Info() << "Reading field p\n" << ref.nl
    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
    
    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
    
    phi = man.createPhi( runTime, mesh, U )
    
    pRefCell = 0
    pRefValue = 0.0
    
    pRefCell, pRefValue = ref.setRefCell( p, mesh.solutionDict().subDict( ref.word( "PIMPLE" ) ), pRefCell, pRefValue )
    
    laminarTransport = man.singlePhaseTransportModel( U, phi )
    
    turbulence = man.incompressible.turbulenceModel.New( U, phi, laminarTransport )
    
    ref.ext_Info() << "Reading field rAU if present\n" << ref.nl
    rAU = man.volScalarField( man.IOobject( ref.word( "rAU" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.READ_IF_PRESENT,
                                            ref.IOobject.AUTO_WRITE ),
                              mesh,
                              runTime.deltaT(),
                              ref.zeroGradientFvPatchScalarField.typeName )
    
    return p, U, phi, laminarTransport, turbulence, rAU, pRefCell, pRefValue


#--------------------------------------------------------------------------------------
def readControls( runTime, mesh, pimple ):
    adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
    
    pimpleDict = pimple.dict()

    correctPhi = pimpleDict.lookupOrDefault( ref.word( "correctPhi" ), ref.Switch( False ) )

    checkMeshCourantNo = pimpleDict.lookupOrDefault( ref.word( "checkMeshCourantNo" ), ref.Switch( False ) )

    return adjustTimeStep, maxCo, maxDeltaT, pimpleDict, correctPhi, checkMeshCourantNo


#--------------------------------------------------------------------------------------
def _correctPhi( runTime, mesh, pimple, p, U, rAU, phi, pRefCell, pRefValue, cumulativeContErr ):
    if mesh.changing():
       for patchi in range( U.ext_boundaryField().size() ):
           if U.ext_boundaryField()[patchi].fixesValue():
              U.ext_boundaryField()[patchi].initEvaluate()
              pass
           pass
       for patchi in range( U.ext_boundaryField().size() ):
           if U.ext_boundaryField()[patchi].fixesValue():
              U.ext_boundaryField()[patchi].evaluate()
              phi.ext_boundaryField()[patchi] << ( U.ext_boundaryField()[patchi] & mesh.Sf().ext_boundaryField()[patchi] )
              pass
           pass
       pass
       
    pcorrTypes = ref.wordList( p.ext_boundaryField().size(), ref.zeroGradientFvPatchScalarField.typeName )
    
    for i in range( p.ext_boundaryField().size() ):
        if p.ext_boundaryField()[i].fixesValue():
           pcorrTypes[i] = ref.fixedValueFvPatchScalarField.typeName
           pass
        pass
    
    pcorr = ref.volScalarField( ref.IOobject( ref.word( "pcorr" ),
                                              ref.fileName( runTime.timeName() ),
                                              mesh,
                                              ref.IOobject.NO_READ,
                                              ref.IOobject.NO_WRITE ),
                                mesh,
                                ref.dimensionedScalar( ref.word( "pcorr" ), p.dimensions(), 0.0),
                                pcorrTypes )
     
    for nonOrth in range( pimple.nNonOrthCorr() + 1 ):
        pcorrEqn = ( ref.fvm.laplacian( rAU, pcorr ) == ref.fvc.div( phi ) )

        pcorrEqn.setReference(pRefCell, pRefValue)
        pcorrEqn.solve()

        if nonOrth == pimple.nNonOrthCorr():
           phi << phi() - pcorrEqn.flux() # mixed calculations
           pass
        pass
    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )     

    return cumulativeContErr

#--------------------------------------------------------------------------------------
def fun_UEqn( mesh, phi, U, p, rAU, turbulence, pimple ):
    # The initial C++ expression does not work properly, because of
    #  1. turbulence.divDevRhoReff( U ) - changes values for the U boundaries
    #  2. the order of expression arguments computation differs with C++
    # UEqn = man.fvm.ddt(U) + man.fvm.div(phi, U) + man.fvVectorMatrix( turbulence.divDevReff( U ), man.Deps( turbulence, U ) ) 
    UEqn = man.fvVectorMatrix( turbulence.divDevReff( U ), man.Deps( turbulence, U ) )  + ( man.fvm.ddt(U) + man.fvm.div(phi, U) )
    
    UEqn.relax()
    
    rAU << 1.0 / UEqn.A()

    if pimple.momentumPredictor():
        ref.solve( UEqn() == -ref.fvc.grad( p ) )
        pass
    else:
        U << rAU *( UEqn.H() - ref.fvc.grad( p ) )
        U.correctBoundaryConditions()
        pass

    return UEqn


#--------------------------------------------------------------------------------------
def  fun_pEqn( mesh, runTime, pimple, U, phi, turbulence, p, rAU, UEqn, pRefCell, pRefValue, cumulativeContErr, corr ):
     U << rAU() * UEqn.H()

     if pimple.nCorr() <= 1 :
        #U Eqn.clear()
        pass
     
     phi << ( ref.fvc.interpolate( U ) & mesh.Sf() )

     if p.needReference():
         ref.fvc.makeRelative( phi, U )
         ref.adjustPhi( phi, U, p )
         ref.fvc.makeAbsolute( phi, U )
         pass
     
     for nonOrth in range( pimple.nNonOrthCorr() + 1 ):
         pEqn = ref.fvm.laplacian( rAU, p ) == ref.fvc.div( phi )
         pEqn.setReference(pRefCell, pRefValue)
         
         pEqn.solve( mesh.solver( p.select( pimple.finalInnerIter( corr, nonOrth ) ) ) )

         if nonOrth == pimple.nNonOrthCorr():
             phi -= pEqn.flux()
             pass
         pass
     cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )

     # Explicitly relax pressure for momentum corrector
     p.relax()
      
     # Make the fluxes relative to the mesh motion
     ref.fvc.makeRelative( phi, U )

     U -= rAU * ref.fvc.grad( p )
     U.correctBoundaryConditions()
     
     return cumulativeContErr

#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createDynamicFvMesh( runTime )
    
    cumulativeContErr = ref.initContinuityErrs()

    p, U, phi, laminarTransport, turbulence, rAU, pRefCell, pRefValue = _createFields( runTime, mesh )    
    
    adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
    
    pimple = man.pimpleControl( mesh )
    
    ref.ext_Info() << "\nStarting time loop\n" <<ref.nl
    
    while runTime.run() :
        
        adjustTimeStep, maxCo, maxDeltaT, pimpleDic, correctPhi, checkMeshCourantNo = readControls( runTime, mesh, pimple )
        CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )
        
        # Make the fluxes absolute
        ref.fvc.makeAbsolute(phi, U)

        runTime = ref.setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )
        runTime.increment()
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        mesh.update()
        if mesh.changing() and correctPhi :
           cumulativeContErr = _correctPhi( runTime, mesh, pimple, p, U, rAU, phi, pRefCell, pRefValue, cumulativeContErr  )
           pass
        
        # Make the fluxes relative to the mesh motion
        ref.fvc.makeRelative( phi, U )
        
        if mesh.changing() and checkMeshCourantNo :
           meshCoNum, meanMeshCoNum = ref.meshCourantNo( runTime, mesh, phi )
           pass
        
        # --- Pressure-velocity PIMPLE corrector loop
        pimple.start()
        while pimple.loop():
            if pimple.nOuterCorr() != 1:
                p.storePrevIter()
                pass
            UEqn = fun_UEqn( mesh, phi, U, p, rAU, turbulence, pimple )
            # --- PISO loop
            for corr in range( pimple.nCorr() ):
                cumulativeContErr = fun_pEqn( mesh, runTime, pimple, U, phi, turbulence, p, rAU, UEqn, pRefCell, pRefValue, cumulativeContErr, corr )
                pass
            
            if pimple.turbCorr():
                turbulence.correct()
                pass
            
            pimple.increment()
            pass
        runTime.write()
        
        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        pass

    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------

from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
else :
   ref.ext_Info() << "\n\n To use this solver it is necessary to SWIG OpenFOAM-2.0.0\n"
   pass

    
#--------------------------------------------------------------------------------------


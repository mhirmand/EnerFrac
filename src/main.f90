
! ---------------------------------------------------------------------
! -                          EnerDyn                                  -
! -       Non-differentiable energy Finite Element Program            -
! -                  for generalized fracture                         -
! - Reference: 
     !	* Hirmand, M. Reza, and Katerina D. Papoulia.              !
     !	 Block coordinate descent energy minimization for 	       !
     !	  dynamic cohesive fracture.   							   !
     !	  Computer Methods in Applied Mechanics and Engineering    !
     !	 354 (2019): 663-688.                                      !
! -                                                                   -
! - Program features:                                                 -
! -   1. Modular structure.                                           -    
! -   2. Explicit time integration with variable step size.           -
! -   3. Hexahedral element with one point gauss qudrature.           -
! -   4. Hourglass control (standard and Flangan-Belytschko methods). -
! -   5. Material model (elasticity, elastic-perfectly plastic,       -
! -      isotropic hardening and Johnson-Cook plasticity models).     -
! -   6. Mie-Gruneisen equation of state.                             -
! -   7. Export results to ParaView for creating and animation and    -
! -      time history curve.                                          -
! -   8. Rigid plane contact                                          -
! -   9. External nodal load and velocity                             -
! -  10. Jaumman rate of stress                                       -
! -                                                                   -
! ---------------------------------------------------------------------

	program EnerDyn
    use constants
	use DataIn
	use DataOut 
	use FFI
    use PreProcess
    !use vtk_generator

	implicit none
	real:: time_begin, time_end
!
	call cpu_time( time_begin )
!	
	call InputPara()
    
    ! update the mesh (duplicate nodes and create interface elements.
    call InitializeData()
!
    !call PreProcess() ! duplicate nodes and create interface elements
!
	call cpu_time( time_end )
	print *, '** Time for preprocessing is', time_end - time_begin, ' seconds'

	call cpu_time( time_begin )
	print *, 'solving...'

	call Integration()	
	
	call cpu_time( time_end )
	print *, '** Time for computation is  ', time_end - time_begin, ' seconds'

	close(ior)
	close(iow1)
	close(iow2)
    close(iow3)

	close(iores)

	end program EnerDyn

	subroutine integration()
! -------------------------------------------------------------------
! -                                                                 -
!   Integrate momentum equations
! -                                                                 -
! -------------------------------------------------------------------
	use Simulation
	use ElementData
	use ParticleData
	use DataOut
	use TimeFunctionData
	implicit none

	real(8) :: prt = 0.0
	real(8) :: plt = 0.0

	real(8) :: tf

	do while (CurrentTime .le. EndTime)
        
		DTk  = DTk1
		DTk1 = 1.0d6

        !	Calculate kinetic energy
		call KineticE()
        
!	Time function at the current time step
		tf = CurrentTMFunction(CurrentTime)	
    
		call FEForceSolids(tf)
        call FEForceFaces(CurrentTime)
            
		call UpdateFEGeometry()
        
        call UpdateDelta(tf)

		if (CurrentTime.ge.prt) then
			prt = prt + OutTimeStep
			call OutCurve(CurrentTime)
			call OutAnim(CurrentTime)
            call write_vtk_file_ascii(istep)
		end if

		if (CurrentTime.ge.plt) then
			plt = plt + RepTimeStep
			write(*,100) istep, CurrentTime, dtk1, EleDistortion, EngKinetic
100			format(1x, "Step=", i6, 1p, "  T=", e10.3, "  DT=", e10.3,  &
					   "  EDist=", e10.3, "  K.E.=", e10.3)
		end if

		istep = istep+1
		CurrentTime = CurrentTime + DTk1

	end do

	end subroutine integration

module MaterialData
    use constants

	type Material
		integer:: MatType		        ! Material type
		integer:: EosType				! EOS type
		real(8):: Density, Young, Poisson
		real(8):: WaveSpeed				! Wave speed c
		real(8):: Yield0, TangMod		! Yield stress and tangent modulus
		real(8):: B_jnc, n_jnc, C_jnc	! constants for Johnson-Cook model
        real(8):: sigma_c, delta_c, G_c ! constants for cohesive fracture model
        
    end type Material
    
    !type CohesiveSurface
    !    integer:: MatType
    !    real(8):: sigma_c, delta_c, G_c ! constants for cohesive fracture model
    !    real(8):: B_jnc, n_jnc, C_jnc	! reserved constants for strain and rate dependent cohesive model
    !end type CohesiveSurface

	integer:: nb_mat = 0				! number of material sets
	type(Material), TARGET, allocatable:: mat_list(:)

	logical:: Jaum  = .true.	! use Jaumman stress rate ?

	type HourGlassControl
		integer:: method = 1	! Type of hourglass control
								!	0 - Off
								!	1 - Standard
								!	2 - Flangan-Belytschko
		real(8):: Qhg = 0.10	! user defined constant for hourglass control
	end type HourGlassControl

	type(HourGlassControl) HourGlass

!	Material history data structure
	type History 
		real(8):: detJ		! Initial volume for finite element
		real(8):: sig_y		! Yield stress
		real(8):: SM, Seqv	! Mean stress and Mises stress
		real(8):: SDS(6)	! deviatoric stress: SDxx, SDyy, SDzz, SDxy, SDyz, SDxz
		real(8):: epeff		! Effective plastic strain
    end type History
    
    !	Material history data structure
	type InterfaceHistory 
		real(8):: sig_max				! maximum traction reached
        real(8):: d_max					! mximum delta reached
        real(8):: jump_u(MAX_DIMENSION) ! jump in u
        real(8):: delta(MAX_DIMENSION)  ! delta
        real(8):: p_avg(3)              ! average traction <sig> n + eta * (jump_u - delta)

	end type InterfaceHistory

	!type(History), TARGET, allocatable:: history_list(:)

! =============================================================================
!   variables used only within the module
!
	private young_, poisson_, yield0_, tangmod_, den0_, epeff_, sig_y_
	private sm, sd, sig, G2, K3, PlaMod, seqv

	real(8):: young_, poisson_	! Young's modulus and Poisson's ratio
	real(8):: yield0_, tangmod_	! Initial Yield stress and tangent modulus
	real(8):: den0_				! Initial density
	real(8):: epeff_, sig_y_	! Effective plastic strain and current yield stress
	real(8):: sm, sd(6), sig(6)	! Mean stress, deviatoric stress and stress components
	real(8):: G2, K3, PlaMod	! 2*G, 3*K, and plastic hardening modulus
	real(8):: seqv				! Equivalent stress

contains


	subroutine InitHistory(historyVars)
! ------------------------------------------------------------------
! -                                                                -
! -	 Purpose                                                       -
! -     initialize history variable for solid point                -
! -                                                                -
!-------------------------------------------------------------------
	implicit none
    type (History), intent(inout):: historyVars
    
	historyVars%sig_y = 0.0d0
	historyVars%SM    = 0.0d0
	historyVars%SDS(1)   = 0.0d0
	historyVars%SDS(2)   = 0.0d0
	historyVars%SDS(3)   = 0.0d0
	historyVars%SDS(4)   = 0.0d0
	historyVars%SDS(5)   = 0.0d0
	historyVars%SDS(6)   = 0.0d0
	historyVars%epeff = 0.0d0

    end subroutine InitHistory

    subroutine initInterfaceHistory(historyVars)
! ------------------------------------------------------------------
! -                                                                -
! -	 Purpose                                                       -
! -     initialize history variables for interface point           -
! -                                                                -
!-------------------------------------------------------------------
	implicit none
    type (InterfaceHistory), intent(inout):: historyVars
    
	historyVars%sig_max = 0.0d0
	historyVars%d_max    = 0.0d0
    historyVars%jump_u(:) = 0.0d0
    historyVars%delta(:)  = 0.0d0
    historyVars%p_avg(:) = 0.0d0
    

	end subroutine initInterfaceHistory

	subroutine Constitution(de, vort, den, his, mat, vol)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -     Material constitution models                                  -
! -        shared by both FE elements and material particles          -
! -                                                                   -
! -  Inputs                                                           -
! -     de		- strain increment                                    -
! -		vort	- vorticity * DT                                      - 
! -     den     - current density                                     -
! -     his  	- history variable                                    -
! -     mat     - material properties                                 -
! -     vol     - current volume                                      -
! -                                                                   -
! -  Return                                                           -
! -		EngStrain - Strain energy                                     -
! -     SDS - Stress components: SDxx, SDyy, SDzz, SDyz, SDxz, SDxy   -
! -     SM, Seqv  - Mean stress and equivalent stress                 -
! -     sig_y     - Yield stress                                      -
! -     epeff     - Effective plastic strain                          -
! -                                                                   -
! -  Note                                                             -
! -		sd(i) and de(i) comply with the Voigt rule                    -
! ---------------------------------------------------------------------
	use Simulation, only: EngStrain, DTk
	implicit none
	real(8), intent(in) :: de(6), vort(3), den, vol
	type(Material), intent(in):: mat 
	type(History), intent(inout):: his

	real(8) dEnS, mu
	integer:: mtype_

!	pick parameters	
	epeff_ = his%epeff
	sig_y_ = his%sig_y

	mtype_ = mat%MatType
	young_ = mat%Young
	den0_  = mat%Density
	poisson_ = mat%Poisson
	yield0_  = mat%Yield0
	tangmod_ = mat%TangMod
	mu = den/den0_ - 1

	G2 = young_ / (1 + poisson_)
	K3 = young_ / (1 - 2*poisson_)
	PlaMod = young_ * tangmod_ / (young_ - tangmod_)	! plastic hardening modulus

!	mean stress
	sm = his%SM

!	deviatoric stress
	sd = his%SDS

!	cauchy stress
	sig(1) = sd(1) + sm
	sig(2) = sd(2) + sm
	sig(3) = sd(3) + sm
	sig(4) = sd(4)
	sig(5) = sd(5)
	sig(6) = sd(6)

!	rotate stress
	if(Jaum) call sigrot(vort, sig, sm, sd)

!	strain energy increment
	dEnS = dot_product(sig, de)

!	select material model
	select case(mtype_)

!	elas: elastic model
    case(1) 
    case(5)
		call M3D1(de)

!	pla1: elastic-perfectly plastic
    case(2) 
    case(6)
		call M3D2(de, den, mu)

!	pla2: isotropic hardening
!	Ref: LS-DYNA Theoretical Manual: Material Model 10
    case(3) 
    case(7)
		call M3D3(de, den, mu)

!	john: Simplfied Johnson-Cook plasticity model
!	Ref: LS-DYNA Theoretical Manual: Material Model 15, 98 (G.R.Johnson 1988)
    case(4) 
    case(8)
		call M3D4(de, den, mu, mat, DTk)

!	Error: invalid material model
	case default
		write(*, 10) mtype_
10		format(1x,'*** Stop *** material type ', i2, ' has not been implemented !')
		stop

	end select

!	Update history variables
	his%SM = sm
	his%Seqv = seqv
	his%SDS = sd
	his%sig_y = sig_y_
	his%epeff = epeff_

!	caculate strain energy and kinetic energy
	dEnS = dEnS + ((sd(1)+sm)*de(1) + (sd(2)+sm)*de(2) + (sd(3)+sm)*de(3)  &
				+ sd(4)*de(4) + sd(5)*de(5) + sd(6)*de(6))

	EngStrain = EngStrain + 0.5d0*dEnS*vol

	end subroutine Constitution


	subroutine sigrot(vort, sig, sm, sd)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Rotate stress if Jaumman stress rate is used                   -
! -     Similar to DYNA3D                                             -
! -                                                                   -
! -  Input                                                            -
! -     vort - Spin * DTk                                             -
! -     sig  - Stress at time Tk                                      -
! -                                                                   -
! -  Output                                                           -
! -     sig  - Updated stress at Tk1                                  -
! -     sm   - Mean stress (pressure)                                 -
! -     sd   - Devitoric stress                                       -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none
	real(8), intent(in) :: vort(3)
	real(8), intent(inout) :: sig(6)
	real(8), intent(out):: sm, sd(6)
	
	real(8) :: rot(6), q(3)
	
	q(1) = 2d0*sig(6)*vort(3)
	q(2) = 2d0*sig(5)*vort(2)
	q(3) = 2d0*sig(4)*vort(1)

	rot(1) =  q(1) + q(2)
	rot(2) = -q(1) + q(3)
	rot(3) = -q(2) - q(3)
	rot(4) = vort(1)*(sig(3)-sig(2)) - vort(2)*sig(6) - vort(3)*sig(5)
	rot(5) = vort(2)*(sig(3)-sig(1)) - vort(1)*sig(6) + vort(3)*sig(4)
	rot(6) = vort(3)*(sig(2)-sig(1)) + vort(1)*sig(5) + vort(2)*sig(4)

	sig = sig + rot

	sm = (sig(1)+sig(2)+sig(3))/3d0

	sd(1) = sig(1) - sm
	sd(2) = sig(2) - sm
	sd(3) = sig(3) - sm
	sd(4) = sig(4)
	sd(5) = sig(5)
	sd(6) = sig(6)
	 
	end subroutine sigrot


	subroutine elastic_devi(de, sd, G2)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Update devitoric stress by elastic relation                   -
! -                                                                   -
! -  Input                                                            -
! -     de - Strain increments                                        -
! -     G2 - 2*G (G is the shear modulus)                             -
! -                                                                   -
! -  Output                                                           -
! -     sd - Devitoric stresses                                       -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none
	
	real(8), intent(in) :: de(6), G2
	real(8), intent(inout) :: sd(6)

	real(8):: dem
	real(8):: dsd(6)

	dem = (de(1) + de(2) + de(3))/3d0
	
	dsd(1) = G2*(de(1)-dem)
	dsd(2) = G2*(de(2)-dem)
	dsd(3) = G2*(de(3)-dem)
	dsd(4) = G2*de(4)*0.5    ! Incremental shear stresses
	dsd(5) = G2*de(5)*0.5
	dsd(6) = G2*de(6)*0.5

	sd = sd + dsd

	end subroutine elastic_devi


	subroutine elastic_p(de, sm, K3)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Update pressure by elastic relation                           -
! -                                                                   -
! -  Input                                                            -
! -     de - Strain increments                                        -
! -     K3 - 3*K (K is the bulk modulus)                              -
! -                                                                   -
! -  Output                                                           -
! -     sm - Mean stress (pressure)                                   -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none
	real(8), intent(in) :: de(6), K3
	real(8), intent(inout) :: sm

	real(8):: dem, dsm

	dem = (de(1) + de(2) + de(3))/3d0
	dsm = K3*dem
	sm = sm + dsm

	end subroutine elastic_p

	real function EquivalentStress(sd)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -     Calculate the equivalent stress sqrt(3*J2)                    -
! -                                                                   -
! -  Input                                                            -
! -     sd	- the deviatoric stress components                        -
! -                                                                   -
! -  Return                                                           -
! -     EquivalentStress - the equivalent stress                      -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	real(8), intent(in) :: sd(6)

	real(8) :: J2	! the second stress invariant

	J2 = 0.5d0*dot_product(sd,sd)
	EquivalentStress = sqrt(J2*3d0)

	return
	end function EquivalentStress


	subroutine M3D1(de)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -     Elastic material model                                        -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	real(8), intent(in) :: de(6)

	call elastic_devi(de, sd, G2)
	call elastic_p(de, sm, K3)
	seqv = EquivalentStress(sd)

	end subroutine M3D1


	subroutine M3D2(de, den, mu)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -     Elastic-perfectly plastic material model                      -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	real(8), intent(in) :: de(6), den, mu

	real(8) ratio, depeff

	call elastic_p(de, sm, K3)

	call elastic_devi(de, sd, G2)
	seqv = EquivalentStress(sd)

	if (seqv .GT. yield0_) then
		ratio = yield0_/seqv
		sd = sd*ratio
		seqv = seqv*ratio
	end if

	end subroutine M3D2


	subroutine M3D3(de, den, mu)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -     Isotropic hardening plastic material model                    -
! -     LS-DYNA theorectical manual : Material Model 10               -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	real(8), intent(in) :: de(6), den, mu

	real(8) ratio, depeff

	call elastic_p(de, sm, K3)

	call elastic_devi(de, sd, G2)
	seqv = EquivalentStress(sd)

	if (seqv .GT. sig_y_) then
		depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
		epeff_ = epeff_ + depeff
		sig_y_ = sig_y_ + PlaMod*depeff

		ratio = sig_y_/seqv
		sd = sd*ratio
		seqv = seqv*ratio
	end if

	end subroutine M3D3


	subroutine M3D4(de, den, mu, mat, DTk)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -     Johnson-Cook plasticity material model (john)                 -
! -     LS-DYNA theorectical manual : Material Model 15               -
! -     G R Johnson 1988                                              -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	real(8), intent(in) :: de(6), den, mu, DTk
	type(Material), intent(in) :: mat

	real(8) ratio, depeff
	real(8):: B_jnc_, n_jnc_, C_jnc_

	call elastic_p(de, sm, K3)

	call elastic_devi(de, sd, G2)
	seqv = EquivalentStress(sd)

	if (seqv .GT. sig_y_) then
		depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
		epeff_ = epeff_ + depeff

		B_jnc_ = mat%B_jnc
		n_jnc_ = mat%n_jnc
		C_jnc_ = mat%C_jnc

		sig_y_ = (yield0_ + B_jnc_*(epeff_**n_jnc_)) * (1 + C_jnc_*log(depeff*1000/DTk))

		ratio = sig_y_/seqv
		sd = sd*ratio
		seqv = seqv*ratio
	end if

	end subroutine M3D4

end module MaterialData
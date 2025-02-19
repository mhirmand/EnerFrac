
	subroutine FEForceSolids(tf)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -    Update FE nodal force                                          -
! -                                                                   -
! -  Input                                                            -
! -    tf - value of time function at current time                    -
! -                                                                   -
! ---------------------------------------------------------------------
	use ParticleData
	use MaterialData
	use ElementData
	use Simulation
    use quadrature
    use constants
	implicit none

	real(8), intent(in) :: tf

	real(8):: vol, den, detJ			! element volume and density 
	real(8):: xyz(3,8), v(3,8)	! element nodal coordinates and velocities 
	real(8):: ijac(3,3)  		! element inverse Jacobian matrix
	real(8):: de(6), vort(3)	! incremental strain and spin
	real(8):: pkxj(4,3)			! derivatives of shape function with respect to coordinates
    real(8), dimension(6, 8*3) :: Bmat
	integer:: ie, igp, i, j, p		! loop counter
    real(8), dimension(8) :: NshapeSolid
    real(8), dimension(8, 3) :: dN_dXi

	real(8) :: dte, at			! element time step and distortion
	real(8) :: c				! element wave speed and hourglass control constant

	type(Element),  POINTER :: el
	type(Material), POINTER :: mat
	type(History),  POINTER :: his

	EleDistortion = 1.0
	
	Acc = Fp*tf
    
    !	loop over all elements
    do ie = 1, nb_element
        el => element_list(ie)

        ! retrieve the nodal coordinates and velocities of the element
        do j = 1, 8
            p = el%nNode(j)
            xyz(:,j) = Pos(:,p)		! Element nodal coordinates
            v(:,j) = Vel(:,p)		! Element nodal velocities
        end do

        ! loop over integration points
        vol = 0.0d0
        do igp = 1, el%numIntPoints

            his => element_list(ie)%historyVars(igp) !history_list(el%nHistory)	! element history variables

            ! evaluate the element inverse Jacobian matrix and volume
            call solidJacobian(xyz, el%numNodes, quadratures_3D(igp,:), ijac, detJ)

            ! the element is severely distorted
            if (detJ .le. 0) then
                write(*, 100) ie
100             format(1x, "*** Error *** Jacobian determinant of element ", i6, " is negative ! ")
                stop
            end if
            
            call solidShapeFunctions(quadratures_3D(igp,:), el%numNodes, NshapeSolid, dN_dXi)

            !	evaluate the increament strain and spin
            call Strain(DTk, el%numNodes, ijac, v, quadratures_3D(igp,:), Bmat, de, vort)

            mat => mat_list(el%mat)				! element material data
            den = mat%density*his%detJ/detJ		! Current density
            c   = mat%WaveSpeed					! wave speed

            !		update stress by constitution law
            call Constitution(de, vort, den, his, mat, detJ)

            ! evaluate and accumulate the equivalent nodal forces due to element stress
            call NodalForceSolid(el%nNode, el%numNodes, Bmat(:, 1:3*el%numNodes), his, detJ * weights_3D(igp))
            
            vol = vol + detJ * weights_3D(igp)

            ! evaluate and accumulate the hourglass-resisting forces
            !if (el%numNodes .eq. 8 .and. el%numIntPoints == 1) then
                ! check the implementation of this subroutine
                !call HourGlassForce(el%nNode, v, xyz, pkxj, vol, den, c)
            !endif
        end do
        
        !		evaluate element time step size and distortion
        call ElementTimeStep(xyz, vol, c, DTk1, EleDistortion)

    end do

    !	determine the new time step size
    if (abs(DTk1-1.0e6) .gt. 1.0e-5) then
        DTk1 = DTk1*DTScale
        DTk1 = min(DTk1, 1.05*DTk1Old)
        DTk1Old = DTk1
    end if

    return
    end subroutine FEForceSolids


	subroutine NodalForce(nNode, pkxj, his, vol)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    Calculate and accumulate the nodal forces due to element stress    -
! -                                                                       -
! - Input                                                                 -
! -    nNode - element nodes                                              -
! -    pkxj  - derivatives of shape function with respect to coordinates  -
! -    his   - element history variable                                   -
! -    vol   - element volume                                             -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
    use constants
	implicit none

	real(8), intent(in) :: pkxj(4,3), vol
	integer, intent(in) :: nNode(8)
	type(History), intent(in) :: his

	real(8) :: sm, sig(6)	! mean stress and stresses 
	real(8) :: f(3,8)		! nodal forces due to element stresses
	integer :: i, k, p		! loop counter

	sm = his%SM  ! Mean stress

	sig(1) = (his%SDS(1) + sm)*vol
	sig(2) = (his%SDS(2) + sm)*vol
	sig(3) = (his%SDS(3) + sm)*vol
	sig(4) = (his%SDS(4))*vol
	sig(5) = (his%SDS(5))*vol
	sig(6) = (his%SDS(6))*vol

!	evaluate nodal forces
	do k = 1, 4
		f(1,k) = -(sig(1)*pkxj(k,1) + sig(6)*pkxj(k,2) + sig(5)*pkxj(k,3))
		f(2,k) = -(sig(2)*pkxj(k,2) + sig(6)*pkxj(k,1) + sig(4)*pkxj(k,3))
		f(3,k) = -(sig(3)*pkxj(k,3) + sig(4)*pkxj(k,2) + sig(5)*pkxj(k,1))
	end do

	do i = 1, 3
		f(i,5) = -f(i,3)
		f(i,6) = -f(i,4)
		f(i,7) = -f(i,1)
		f(i,8) = -f(i,2)
	end do

!	accumulate nodal forces
	do k = 1, 8
		p = nNode(k)
		Acc(:,p) = Acc(:,p) + f(:,k)
	end do

    end subroutine NodalForce
    
    subroutine NodalForceSolid(nNode, numNodes, Bmat, his, vol)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    Calculate and accumulate the nodal forces due to element stress    -
! -                                                                       -
! - Input                                                                 -
! -    nNode - nodes of the element                                       -
! -    BmatI - B matrix                                                   -
! -    his   - element history variable                                   -
! -    vol   - element volume                                             -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
    use constants
	implicit none

	real(8), intent(in) :: Bmat(6,3*numNodes), vol !
	integer, intent(in) :: numNodes, nNode(numNodes)
	type(History), intent(in) :: his

	real(8) :: sm, sig(6)	! mean stress and stresses 
	real(8) :: fI(3*numNodes)		! nodal forces due to element stresses at node I
	integer :: i, k, p		! loop counter

	sm = his%SM  ! Mean stress

	sig(1) = (his%SDS(1) + sm)*vol
	sig(2) = (his%SDS(2) + sm)*vol
	sig(3) = (his%SDS(3) + sm)*vol
	sig(4) = (his%SDS(4))*vol
	sig(5) = (his%SDS(5))*vol
	sig(6) = (his%SDS(6))*vol

!	evaluate nodal forces
    fI = matmul(transpose(Bmat), sig)
	!do k = 1, 6
	!	fI(1) = fI(1) - BmatI(k, 1) * sig(k) !-(sig(1)*pkxj(k,1) + sig(6)*pkxj(k,2) + sig(5)*pkxj(k,3))
	!	fI(2) = fI(2) - BmatI(k, 2) * sig(k)
	!	fI(3) = fI(3) - BmatI(k, 3) * sig(k)
	!end do

!	accumulate nodal forces
	do k = 1, numNodes
		p = nNode(k)
        Acc(:, p) = Acc(:, p) - fI(3*k-2:3*k)
	end do

    end subroutine NodalForceSolid
    
    subroutine NodalForceFace(nNode, numNodes, NshapeSolid, traction, area)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    Calculate and accumulate the nodal forces due to element stress    -
!      on its face
! -                                                                       -
! - Input                                                                 -
! -    nNode   - face nodes                                               -
! -    Nshape  - shape function on the face                               -
! -    his     - face history variable                                    -
! -    area    - area of the GP (detJ * weight)                           -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
    use constants
	implicit none

    integer, intent(in) :: numNodes
    integer, intent(in) :: nNode(numNodes)
	real(8), intent(in) :: NshapeSolid(numNodes), area


    real(8), intent(in) :: traction(3)	! traction on the face
	real(8) :: f(3,8)		! nodal forces due to element stresses
	integer :: i, k, p		! loop counter


!	evaluate nodal forces
    f(:,:) = 0.0
	do k = 1, numNodes
		f(1,k) = NshapeSolid(k) * traction(1) * area
		f(2,k) = NshapeSolid(k) * traction(2) * area
        f(3,k) = NshapeSolid(k) * traction(3) * area
    end do
    
!	accumulate nodal forces
	do k = 1, numNodes
		p = nNode(k)
		Acc(:,p) = Acc(:,p) - f(:,k)
	end do

	end subroutine NodalForceFace



	subroutine UpdateFEGeometry()
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    Compute accelerations, apply displacement b.c.'s, and update       -
! -       velocity and position for FE nodes                              -
! -                                                                       -
! -------------------------------------------------------------------------
	use Simulation
	use ElementData
	use ParticleData
	implicit none

	integer p, nDir, nUnit
	real(8) :: DT2

	DT2 = (Dtk+DTk1)*0.5

!	Compute accelerations of FE nodes

	do p = 1, nb_node
		Acc(:,p) = Acc(:,p) / Mp(p)
	end do

!	Apply displacement boundary conditions - rigid wall

	nDir = plane%nDir		! unit normal vector of the rigid wall

	if (nDir .ne. 0) then
		nUnit = sign(1, nDir)

		do p = 1, nb_node

!			the rigid wall is perpendicular to x axis
			if (abs(nDir).eq.1 .and. nUnit*(plane%coor-Pos(1,p)).ge.0) then
				if (nUnit*Vel(1,p).lt.0) Vel(1,p) = 0
				if (nUnit*Acc(1,p).lt.0) Acc(1,p)  = 0
			end if

!			the rigid wall is perpendicular to y axis
			if (abs(nDir).eq.2 .and. nUnit*(plane%coor-Pos(2,p)).ge.0) then
				if (nUnit*Vel(2,p).lt.0) Vel(2,p) = 0
				if (nUnit*Acc(2,p).lt.0) Acc(2,p)  = 0
			end if

!			the rigid wall is perpendicular to z axis
			if (abs(nDir).eq.3 .and. nUnit*(plane%coor-Pos(3,p)).ge.0) then
				if (nUnit*Vel(3,p).lt.0) Vel(3,p) = 0
				if (nUnit*Acc(3,p).lt.0) Acc(3,p)  = 0
			end if

		end do
	end if

!	Update velocity and position for FE nodes
	Vel = Vel + Acc*DT2

    end subroutine UpdateFEGeometry
    
    subroutine applyVelocityBCs()
    ! -------------------------------------------------------------------------
    ! - Purpose                                                               -
    ! -    Applies velocity boundary conditions                               -
    ! -                                                                       -
    ! -------------------------------------------------------------------------
	use Simulation
	use ElementData
	use ParticleData
    use constants
	implicit none
    
    integer:: i
    
    
    do i = 1, nb_node
        if (VelBCs(i)%directionIDs(1) == 1) then
            Vel(1,i) = VelBCs(i)%bcVals(1)
            Acc(1,i) = 0.0
        endif
        if (VelBCs(i)%directionIDs(2) == 1) then
            Vel(2,i) = VelBCs(i)%bcVals(2)
            Acc(2,i) = 0.0
        end if
        if (VelBCs(i)%directionIDs(3) == 1) then
            Vel(3,i) = VelBCs(i)%bcVals(3)
            Acc(3,i) = 0.0
        end if
    end do
    end subroutine applyVelocityBCs

	subroutine HourGlassForce(nNode, v, xyz, pkxj, vol, den, c)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate and accumulate hourglass-resisting nodal forces           -
! -                                                                       -
! - Input                                                                 -
! -    nNode - element nodes                                              -
! -    v     - element nodeal velocities                                  -
! -    xyz   - same as before                                             -
! -    pkxj  - same as before                                             -
! -    vol   - element volume                                             -
! -    den   - element density                                            -
! -    c     - material sound speed                                       -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
    use constants
	implicit none

	integer, intent(in) :: nNode(8)
	real(8), intent(in) :: v(3,8), xyz(3,8), pkxj(4,3) 
	real(8), intent(in) :: vol, den, c

	real(8) :: v3478(3), v2358(3), v1467(3), v1256(3)
	real(8) :: hgr(3,4)		! projection of velocity filed on to hourglass modes 
	real(8) :: f(3,8)		! coefficient and hourglass-resisting forces

	real(8) :: hap(3), ham(3), hbp(3), hbm(3)

	real(8) :: a, v32, ah, al
	integer :: i, k, p		! loop counter

	if (HourGlass%method.eq.0 .or. HourGlass%Qhg.le.1.0e-6) return

	if (HourGlass%method .eq. 1) then
		call HGForceStandard(v, vol, den, c, f)
	else if (HourGlass%method .eq. 2) then
		call HGForceFlangan(v, xyz, pkxj, vol, den, c, f)
	else
		stop "***Error*** Invalid Hourglass control method !"
	end if

!	accumulate the hourglass-resisting forces
	do k = 1, 8
		p = nNode(k)
		Acc(:,p) = Acc(:,p) + f(:,k)
	end do

	end subroutine HourGlassForce


	subroutine HGForceStandard(v, vol, den, c, f)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate hourglass-resisting nodal forces using standard method    -
! -                                                                       -
! - Input                                                                 -
! -    v     - element nodeal velocities                                  -
! -    vol   - element volume                                             -
! -    den   - element density                                            -
! -    c     - material sound speed                                       -
! -                                                                       -
! - Output                                                                -
! -    f     - hourglass-resisting nodal forces                           -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
    use constants
	implicit none

	real(8), intent(in) :: v(3,8), vol, den, c
	real(8), intent(out):: f(3,8)

	real(8) :: v3478(3), v2358(3), v1467(3), v1256(3)
	real(8) :: hgr(3,4)		! projection of velocity filed on to hourglass modes 

	real(8) :: hap(3), ham(3), hbp(3), hbm(3)

	real(8) :: qh, a, v32, ah, al
	integer :: i, k, p		! loop counter

	qh = HourGlass%Qhg

	a  = qh*den/4.0
	v32= vol**(2.0/3.0)
	ah = a*v32*c
	al = a*v32*(100*qh)

	do i = 1, 3
		v3478(i) = v(i,3) - v(i,4) - v(i,7) + v(i,8)
		v2358(i) = v(i,2) - v(i,3) - v(i,5) + v(i,8)
		v1467(i) = v(i,1) - v(i,4) - v(i,6) + v(i,7)
		v1256(i) = v(i,1) - v(i,2) - v(i,5) + v(i,6)
	end do

	do i = 1, 3
		hgr(i,1) = v1467(i) - v2358(i)
		hgr(i,2) = v1467(i) + v2358(i)
		hgr(i,3) = v1256(i) - v3478(i)
		hgr(i,4) = v1256(i) + v3478(i)
	end do

	do i = 1, 3
		do k = 1, 4
			hgr(i,k) = hgr(i,k)*(ah + abs(al*hgr(i,k)))
		end do
	end do

	do i = 1, 3
		hap(i) = hgr(i,1) + hgr(i,4)
		ham(i) = hgr(i,1) - hgr(i,4)
		hbp(i) = hgr(i,2) + hgr(i,3)
		hbm(i) = hgr(i,2) - hgr(i,3)
	end do

!	evaluate the hourglass-resisting forces
	do i = 1, 3
		f(i,1) = -hap(i) - hbp(i)
		f(i,2) =  hap(i) - hbm(i)
		f(i,3) = -hap(i) + hbp(i)
		f(i,4) =  hap(i) + hbm(i)
		f(i,5) = -ham(i) + hbp(i)
		f(i,6) =  ham(i) + hbm(i)
		f(i,7) = -ham(i) - hbp(i)
		f(i,8) =  ham(i) - hbm(i)
	end do

	end subroutine HGForceStandard


	subroutine HGForceFlangan(v, xyz, pkxj, vol, den, c, f)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate hourglass-resisting nodal forces using Flangan-Belytschko -
! -       method                                                          -
! -                                                                       -
! - Input                                                                 -
! -    v     - element nodeal velocities                                  -
! -    vol   - element volume                                             -
! -    den   - element density                                            -
! -    c     - material sound speed                                       -
! -                                                                       -
! - Output                                                                -
! -    f     - hourglass-resisting nodal forces                           -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
    use constants
	implicit none

	real(8), intent(in) :: v(3,8), xyz(3,8), pkxj(4,3) 
	real(8), intent(in) :: vol, den, c
	real(8), intent(out):: f(3,8)

	integer, parameter :: h(4,4) = (/ 1,-1, 1,-1, & 
									  1, 1,-1,-1, &
									  1,-1,-1, 1, &
									  1,-1, 1,-1  /)

	real(8), parameter :: ss(4,4) = (/ 2.0,-2.0, 2.0,-2.0,  &
									  -2.0,-2.0, 2.0, 2.0,  &
									  -2.0, 2.0, 2.0,-2.0,  &
									   0.0, 0.0, 0.0, 0.0  /)

	real(8) :: x3478(3), x2358(3), x1467(3), x1256(3)
	real(8) :: hgr(3,4), beta(3,4), gama(4,8), g(3,4)

	real(8) :: ah, qh
	integer :: i, j, k, p		! loop counter

	qh = HourGlass%Qhg

	ah = -qh*c*den*vol**(2.0/3.0)/4.0

	do i = 1, 3
		x3478(i) = xyz(i,3) - xyz(i,4) - xyz(i,7) + xyz(i,8)
		x2358(i) = xyz(i,2) - xyz(i,3) - xyz(i,5) + xyz(i,8)
		x1467(i) = xyz(i,1) - xyz(i,4) - xyz(i,6) + xyz(i,7)
		x1256(i) = xyz(i,1) - xyz(i,2) - xyz(i,5) + xyz(i,6)
	end do

	do i = 1, 3
		beta(i,1) = x1467(i) - x2358(i)
		beta(i,2) = x1467(i) + x2358(i)
		beta(i,3) = x1256(i) - x3478(i)
		beta(i,4) = x1256(i) + x3478(i)
	end do

	do j = 1, 4
		do k = 1, 4
			gama(j,k) = h(k,j)
			do i = 1, 3
				gama(j,k) = gama(j,k) - beta(i,j)*pkxj(k,i) 
			end do
		end do 
	end do

	do j = 1, 4
		gama(j,5) = ss(1,j) - gama(j,3)
		gama(j,6) = ss(2,j) - gama(j,4)
		gama(j,7) = ss(3,j) - gama(j,1)
		gama(j,8) = ss(4,j) - gama(j,2)
	end do

	do i = 1, 3
		do j= 1, 4
			g(i,j) = 0.0
			do k = 1, 8
				g(i,j) = g(i,j) + v(i,k)*gama(j,k)
			end do
		end do
	end do

!	evaluate the hourglass-resisting forces
	do i = 1, 3
		do k = 1, 8
			f(i,k) = 0.0
			do j = 1, 4
				f(i,k) = f(i,k) + g(i,j)*gama(j,k)
			end do
			f(i,k) = ah*f(i,k)
		end do
	end do

	end subroutine HGForceFlangan


	subroutine solidJacobian_old(xyz, vol, ijac)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -    Evaluate Jacobian determinant and inverse Jacobian matrix at   -
! -       the center of the element                                   -
! -                                                                   -
! -  Input                                                            -
! -    xyz  - Coordinates of vertexes of the element                  -
! -                                                                   -
! -  Output                                                           -
! -    vol  - Volume of the element                                   -
! -    ijac - Inverse Jacobian matrix                                 -
! -                                                                   -
! ---------------------------------------------------------------------
!	use lin_sol_gen_int		! IMSL Fortran 90 MP LIBRARY 
							! need to link with SF90MP.LIB SMATHD.LIB
	!use imsl
    use constants
	implicit none
    

	real(8), intent(in)  :: xyz(3,8)  !  Coordinates of vertexes of the element
	real(8), intent(out) :: vol       !  Volume of the element
	real(8), intent(out) :: ijac(3,3) !  Inverse Jacobian

	real(8):: det_       !  Determinant of jacobian
	real(8):: jac(3,3)  !  Jacobian matrix

	real(8):: x17(3), x28(3), x35(3), x46(3)
	real(8):: b(3,0), x(3,0), dett(2)   !  Required for lin_sol_gen
	integer:: j

!   JACOBIAN MATRIX at the center of a element (one-point quadrature)
	do j = 1, 3
		x17(j) = xyz(j,7) - xyz(j,1)
		x28(j) = xyz(j,8) - xyz(j,2)
		x35(j) = xyz(j,5) - xyz(j,3)
		x46(j) = xyz(j,6) - xyz(j,4)
	end do

	do j = 1, 3
		jac(1,j) = (x17(j) - x28(j) - x35(j) + x46(j))/8.0
		jac(2,j) = (x17(j) + x28(j) - x35(j) - x46(j))/8.0
		jac(3,j) = (x17(j) + x28(j) + x35(j) + x46(j))/8.0
	enddo

!	Compute the matrix inverse and its determinant. 
 	!call lin_sol_gen(jac, b, x, nrhs=0, ainv=ijac, det=dett) 

	!det_ = abs(dett(1))**dett(2) * (dett(1))/abs(dett(1))
    
    call innve(jac,ijac,det_)

!	element volume
	vol = 8.0*det_

	return
    end subroutine solidJacobian_old
    
    subroutine solidJacobian(nodalCoords, numNodes, natCoord, invJacobian, detJ)
    ! ---------------------------------------------------------------------
	! -  Purpose                                                          -
	! -    Evaluate Jacobian determinant and inverse Jacobian matrix at   -
	! -       a given point in the element                                -
	! -                                                                   -
	! -  Input                                                            -
	! -    xyz  - Coordinates of vertexes of the element                  -
	! -                                                                   -
	! -  Output                                                           -
	! -    vol  - Volume of the element                                   -
	! -    ijac - Inverse Jacobian matrix                                 -
	! -                                                                   -
! ---------------------------------------------------------------------
  use constants

  implicit none
  real(8), dimension(3, numNodes), intent(in) :: nodalCoords  ! Coordinates of the 8 nodes (x, y, z)
  real(8), dimension(3), intent(in)   :: natCoord     ! Natural coordinates (xi, eta, zeta)
  real(8), dimension(3,3), intent(out):: invJacobian  ! inverse of Jacobian matrix
  real(8), intent(out)                :: detJ         ! Determinant of the Jacobian
  integer, intent(in)                 :: numNodes     ! Number of nodes in the element

  real(8), dimension(numNodes,3) :: dN_dXi             ! Shape function derivatives w.r.t. natural coordinates
  real(8), dimension(numNodes) :: shapeFunc                    ! Shape functions
  real(8) :: jacobian(3,3)                                   !  Jacobian matrix
  integer :: i, j
  real(8) :: detJ2

  ! Initialize Jacobian to zero
  jacobian = 0.0d0
  invJacobian = 0.0d0
  
  ! evaluate the shape function derivatives with respect to natural coordinates
  call solidShapeFunctions(natCoord, numNodes, shapeFunc, dN_dXi)

  ! Compute the Jacobian matrix
  do i = 1, 3
    do j = 1, 3
      jacobian(i,j) = sum(dN_dXi(:,j) * nodalCoords(i, :))
    end do
  end do

  ! Compute the determinant of the Jacobian matrix
  detJ2 = jacobian(1,1) * (jacobian(2,2) * jacobian(3,3) - jacobian(2,3) * jacobian(3,2)) &
       - jacobian(1,2) * (jacobian(2,1) * jacobian(3,3) - jacobian(2,3) * jacobian(3,1)) &
       + jacobian(1,3) * (jacobian(2,1) * jacobian(3,2) - jacobian(2,2) * jacobian(3,1))
  
  ! compute the inverse of the Jacobian matrix
  call innve(jacobian,invJacobian, detJ)
  
  
end subroutine solidJacobian

subroutine Strain_old(DTk, v, ijac, de, vort, pkxj)
! ---------------------------------------------------------------------
! - Purpose                                                           -
! -    Calculate incremental strain and vorts                         -
! -                                                                   -
! - Input                                                             -
! -    DTk   - Time step size                                         -
! -    v     - element nodal velocities                               -
! -    ijac  - Inverse Jacobian                                       -
! -                                                                   -
! - Output                                                            -
! -    de    - Velocity strain rates                                  -
! -    vort  - Non zero entries of incremental spin tensor            -
! -    pkxj  - derivatives of shape function d(Phi_k)/d(x_j)          -
! -                                                                   -
! ---------------------------------------------------------------------
    use constants	
    implicit none

	real(8), intent(in)  :: DTk        ! Current time step size
	real(8), intent(in)  :: v(3,8)     ! element nodal velocity
	real(8), intent(in)  :: ijac(3,3)  ! Inverse Jacobian
	real(8), intent(out) :: de(6)      ! Velocity strain rates
	real(8), intent(out) :: vort(3)    ! incremental spin
	real(8), intent(out) :: pkxj(4,3)  ! derivatives of shape function d(Phi_k)/d(x_j)

	real(8) :: v17(3), v28(3), v35(3), v46(3)
	real(8) :: dij(3,3)
	integer i, j

!   evaluate the derivatives of shape function with respect to spatial coordinates x_j
	do j = 1, 3
		pkxj(1,j) = (-ijac(j,1) - ijac(j,2) - ijac(j,3))/8.0d0
		pkxj(2,j) = ( ijac(j,1) - ijac(j,2) - ijac(j,3))/8.0d0
		pkxj(3,j) = ( ijac(j,1) + ijac(j,2) - ijac(j,3))/8.0d0
		pkxj(4,j) = (-ijac(j,1) + ijac(j,2) - ijac(j,3))/8.0d0
	end do
!
	do i = 1, 3
		v17(i) = v(i,1) - v(i,7)
		v28(i) = v(i,2) - v(i,8)
		v35(i) = v(i,3) - v(i,5)
		v46(i) = v(i,4) - v(i,6)
	end do

	do i = 1, 3
		do j = 1, 3
			dij(i,j) = v17(i)*pkxj(1,j)+v28(i)*pkxj(2,j)+v35(i)*pkxj(3,j)+v46(i)*pkxj(4,j)
		end do
		de(i) = dij(i,i)   ! Normal velocity strain rates
	end do

	de(4) = dij(2,3) + dij(3,2)   ! Engineering shear velocity strain rates
	de(5) = dij(1,3) + dij(3,1)
	de(6) = dij(1,2) + dij(2,1)

	vort(1) = (dij(2,3) - dij(3,2))*0.5   ! spin : wyz
	vort(2) = (dij(1,3) - dij(3,1))*0.5   ! spin : wxz
	vort(3) = (dij(1,2) - dij(2,1))*0.5   ! spin : wxy

	de = de * DTk		! Incremental strain
	vort = vort *DTk	! Incremental vort

end subroutine Strain_old
    
subroutine Strain(dt, numNode, invJacobian, nodalVelocities, intPoint, B, de, vort)
! ---------------------------------------------------------------------
! - Purpose                                                           -
! -    Calculate incremental strain and vorts                         -
! -                                                                   -
! - Input                                                             -
! -    dt                  - Time step size                           -
! -    nodalVelocities     - element nodal velocities                 -
! -    invJacobian         - Inverse Jacobian                         -
! -    intPoint            - natural coords of integration point      -
! -                                                                   -
! - Output                                                            -
! -    de    - Velocity strain rates                                  -
! -    vort  - Non zero entries of incremental spin tensor            -
! -    B     - strain-displacement matrix                             -
! -                                                                   -
! ---------------------------------------------------------------------
  implicit none
  integer, intent(in) :: numNode							! Number of nodes in the element
  real(8), intent(in) :: dt									! Time increment
  real(8), dimension(3,3), intent(in) :: invJacobian		! Inverse Jacobian matrix
  real(8), dimension(3, numNode), intent(in) :: nodalVelocities	! Nodal velocities (vx, vy, vz), each node in one column
  real(8), dimension(3), intent(in) :: intPoint				! Natural coordinates of integration point
  
  real(8), dimension(6, 3 * numNode), intent(out) :: B      ! B matrix (6 rows for strain components, 8 columns for nodes)
  real(8), dimension(6), intent(out) :: de                  ! Strain increment
  real(8), dimension(3), intent(out) :: vort				! Spin increment (vorticity)

  real(8), dimension(3,3) :: L								! Velocity gradient tensor
  real(8), dimension(3,3) :: strainIncrement				! Strain increment tensor
  real(8), dimension(3,3) :: spinIncrement					! Spin increment tensor
  real(8), dimension(numNode,3) :: dN_dX					! Shape function derivatives w.r.t. global coordinates
  real(8), dimension(numNode,3) :: dN_dXi				    ! Shape function derivatives w.r.t. natural coordinates
  real(8), dimension(numNode) :: Nshape							! Shape functions


  integer :: i, j, k

  ! Initialize tensors
  L = 0.0d0
  strainIncrement = 0.0d0
  spinIncrement = 0.0d0
  dN_dX = 0.0d0
  B = 0.0d0
  
  call solidShapeFunctions(intPoint, numNode, Nshape, dN_dXi)
  
! Compute shape function derivatives w.r.t. global coordinates: dN_dX = J_inv * dN_dXi
  do i = 1, 3
    do k = 1, 8
      do j = 1, 3
        dN_dX(k,i) = dN_dX(k,i) + invJacobian(i,j) * dN_dXi(k,j)
      end do
    end do
  end do

  ! Compute velocity gradient tensor L = J_inv * dN_dXi * nodalVelocities
  do i = 1, 3
    do j = 1, 3
      do k = 1, 8
        L(i,j) = L(i,j) + dot_product(invJacobian(i,:), dN_dXi(k,:)) * nodalVelocities(j,k)
      end do
    end do
  end do

  ! Compute strain and spin increments
  do i = 1, 3
    do j = 1, 3
      strainIncrement(i,j) = 0.5d0 * (L(i,j) + L(j,i)) * dt
      spinIncrement(i,j)   = 0.5d0 * (L(i,j) - L(j,i)) * dt
    end do
  end do

  ! put strain and spin increments into vectors 
  de(1:6) = (/ strainIncrement(1,1), strainIncrement(2,2), strainIncrement(3,3), &
                      2.0d0*strainIncrement(1,2), 2.0d0*strainIncrement(1,3), 2.0d0*strainIncrement(2,3) /) ! engineering shear definitions

  vort(1:3) = (/ spinIncrement(1,2), spinIncrement(1,3), spinIncrement(2,3) /)
  
    ! Assemble B matrix (strain-displacement matrix)
  do i = 1, numNode
    B(1, 3*i - 2) = dN_dX(i,1)  ! Strain xx
    B(2, 3*i - 1) = dN_dX(i,2)  ! Strain yy
    B(3, 3*i) = dN_dX(i,3)  ! Strain zz
    !
    B(4, 3*i - 2) = dN_dX(i,2)  ! Strain xy (symmetry)
    B(4, 3*i - 1) = dN_dX(i,1)
    !
    B(5, 3*i - 2) = dN_dX(i,3)  ! Strain xz (symmetry)
    B(5, 3*i) = dN_dX(i,1)
    !
    B(6, 3*i - 1) = dN_dX(i,3)  ! Strain yz (symmetry)
    B(6, 3*i) = dN_dX(i,2)
  end do
  
end subroutine Strain



	subroutine ElementTimeStep(xyz, vol, c, DTk1, EleDistortion)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate element time step size and distortion                     -
! -                                                                       -
! - Input                                                                 -
! -    xyz  - element nodal coordinates                                   -
! -    vol  - element volume                                              -
! -    c    - wave speed                                                  -
! -    DTk1          - current time step size                             -
! -    EleDistortion - current element distortion                         -
! -                                                                       -
! - Output                                                                -
! -    DTk1          - current time step size                             -
! -    EleDistortion - current element distortion                         -
! -                                                                       -
! -------------------------------------------------------------------------
	use constants	
    implicit none

	real(8), intent(in) :: xyz(3,8), vol, c
	real(8), intent(inout) :: DTk1, EleDistortion

	integer :: fac(4,6) = (/ 1,2,3,4, 5,6,7,8, 1,2,6,5, &	! element surface definition
							 2,3,7,6, 3,4,8,7, 4,1,5,8 /)

	real(8) ::dt, at
	real(8) :: e, f, g, atest, areal, aream, x13(3), x24(3), fs(3), ft(3)
	integer :: i, j, k1, k2, k3, k4

	areal = 1.0e20
	aream = 0.0

!	Loop over all six surfaces of the element
	do j = 1, 6
		do i = 1, 3
			k1 = fac(1,j)
			k2 = fac(2,j)
			k3 = fac(3,j)
			k4 = fac(4,j)

			x13(i) = xyz(i,k3) - xyz(i,k1)
			x24(i) = xyz(i,k4) - xyz(i,k2)

			fs(i) = x13(i) - x24(i)
			ft(i) = x13(i) + x24(i)
		end do

		e = fs(1)*fs(1) + fs(2)*fs(2) + fs(3)*fs(3)
		f = fs(1)*ft(1) + fs(2)*ft(2) + fs(3)*ft(3)
		g = ft(1)*ft(1) + ft(2)*ft(2) + ft(3)*ft(3)

		atest = e*g - f*f	! area/4

		aream = max(atest, aream)
		areal = min(atest, areal)

	end do

	at = areal/aream
	dt = 4*vol/sqrt(aream)
	dt = dt/c / sqrt(betaDG)

	DTk1 = min(DTk1, dt)
	EleDistortion = min(at, EleDistortion)

	end subroutine ElementTimeStep


	subroutine KineticE
! ------------------------------------------------------------------
! -                                                                -
! -  purpose:                                                      -
! -     calculate kinematic energy                                 -
! -                                                                -
!-------------------------------------------------------------------
	use ParticleData
	use Simulation
    use constants	
	implicit none

	integer :: p	! loop counter

	EngKinetic = 0.0

	do p = 1, nb_node
		EngKinetic = EngKinetic + dot_product(Vel(:,p), Vel(:,p))*Mp(p)*0.5d0
	end do
    end subroutine KineticE
    
    
    
    subroutine FEForceFaces(tf)
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -    ! Evaluate and accumulate forces due to interfaces             -
    ! -                                                                   -
    ! -  Input                                                            -
    ! -    tf - value of time function at current time                    -
    ! -                                                                   -
    ! ---------------------------------------------------------------------
    use ParticleData
    use MaterialData
    use ElementData
    use Simulation
    use quadrature
    use constants	
    
    implicit none
    real(8), intent(in) :: tf
    integer :: ie, iface, igp, j, p, f, n, m, k
    integer :: numFaces, directionID
    
    real(8):: den1, den2			    ! density 
    integer :: solid1, solid2, face1, face2
    real(8):: xyzFace(3,4), xyzSolid1(3, 8), xyzSolid2(3, 8), v1(3,8), v2(3, 8)	! element nodal coordinates and velocities 
	real(8):: invJacobian1(3,3), invJacobian2(3,3)		! element inverse Jacobian matrix
	real(8):: de1(6), de2(6), vort1(3), vort2(3)	! incremental strain and spin
	! real(8):: pkxj(4,3)			! derivatives of shape function with respect to coordinates
    real(8):: faceDetJ, solidDetJ1,  solidDetJ2            ! determinant of the Jacobian matricies
    real(8), dimension(4, 3) :: quadratures_faces_3D ! quadrature points on the face in 3D natural coords of the element
    real(8), dimension(6, 24) :: B1, B2
    real(8), dimension(4) :: NshapeFace
    real(8), dimension(8) :: NshapeSolid1, NshapeSolid2
    real(8), dimension(3) :: jump_u, delta
    real(8), dimension(8, 3) :: dN_dXi
    real(8), dimension(3) :: faceNormal, traction
    real(8) :: lface, sm
    real(8), dimension(6) :: sig
    

    type(InterfaceElem),  POINTER :: intfel
    type(Element),  POINTER :: solel1, solel2
    type(Face), POINTER :: eface
	type(Material), POINTER :: mat
	type(History),  POINTER :: his1,his2

    ! loop over all interface elements 
    do ie = 1, nb_interface

        ! get intf element
        intfel => interface_list(ie)

        
        ! numFaces = el%numFaces
        
        ! get solid IDS
        solid1 = intfel%solidID1
        solid2 = intfel%solidID2
        
        ! get face ID of the solids
        face1 = intfel%faceID1
        face2 = intfel%faceID2
        
        ! get the solid elements
        solel1 => element_list(solid1)
        solel2 => element_list(solid2)
        
        ! get solid elements coords and velocities
        v1 = 0.0d0
        do j = 1, solel1%numNodes
            p = solel1%nNode(j)
            xyzSolid1(:,j) = Pos(:,p)		! Element nodal coordinates
            v1(:,j) = Vel(:,p)		        ! Element nodal velocities
        end do

        v2 = 0.0d0
        do j = 1, solel2%numNodes
            p = solel2%nNode(j)
            xyzSolid2(:,j) = Pos(:,p)		! Element nodal coordinates
            v2(:,j) = Vel(:,p)		        ! Element nodal velocities
        end do

        ! loop over all faces
        do igp = 1, intfel%numIntPoints
            
            ! delta = 0.0d0
            delta = intfel%historyVars(igp)%delta
            
            ! skip if delta_max > delta_c
            
            ! get the face nodes of the first edge
            do j = 1, solel1%faces(face1)%numNodes
                p = solel1%faces(face1)%nNode(j)
                xyzFace(:,j) = Pos(:,p)		! Element nodal coordinates
            end do

            ! evluate face length
            lface = sqrt(sum((xyzFace(:, 1) - xyzFace(:, 3))**2))

                ! evaluate the element inverse Jacobian matrix on the face for strain evluation
                call solidJacobian(xyzSolid, el%numNodes, eface%intPoints(igp, :), invJacobian, solidDetJ)

            !do igp = 1, eface%numIntPoints
            his1 => solel1%faces(face1)%historyVars(igp) !history_list(el%nHistory)	! element history variables
            his2 => solel2%faces(face2)%historyVars(igp)
                
                ! the element is severely distorted
            if (faceDetJ .le. 0.0d0) then
                write(*, 100) face1, ie
                100 format(1x, "*** Error *** Interface Jacobian determinant of face ", i6, " in element " , i6, " is negative ! ")
                    stop
                end if

            ! evaluate the element inverse Jacobian matrix on the face for strain evluation
            call solidJacobian(xyzSolid1, solel1%numNodes, solel1%faces(face1)%intPoints(igp, :), invJacobian1, solidDetJ1)
            call solidJacobian(xyzSolid2, solel2%numNodes, solel2%faces(face2)%intPoints(igp, :), invJacobian2, solidDetJ2)

            ! the element is severely distorted
            if (SolidDetJ1 .le. 0.0d0) then
                write(*, 101) face1, solid1
                101 format(1x, "*** Error *** Solid Jacobian determinant of face ", i6, " in element " , i6, " is negative !")
                !stop
            end if

            ! the element is severely distorted
            if (SolidDetJ2 .le. 0.0d0) then
                write(*, 102) face2, solid2
                102 format(1x, "*** Error *** Solid Jacobian determinant of face ", i6, " in element " , i6, " is negative !")
                !stop
            end if
                
                ! evaluate the increament strain and spin
                !call Strain_old(DTk, v, invJacobian, de, vort, pkxj)
                call solidShapeFunctions(solel1%faces(face1)%intPoints(igp, :), solel1%numNodes, NshapeSolid1, dN_dXi)
                call solidShapeFunctions(solel2%faces(face2)%intPoints(igp, :), solel2%numNodes, NshapeSolid2, dN_dXi)
                
                ! evaluate the jump in displacement
                jump_u = matmul(xyzSolid2, NshapeSolid2) - matmul(xyzSolid1, NshapeSolid1)
                
                !evaluate strain
                call Strain(DTk, solel1%numNodes, invJacobian1, v1, solel1%faces(face1)%intPoints(igp, :), B1, de1, vort1)
                call Strain(DTk, solel2%numNodes, invJacobian2, v2, solel2%faces(face2)%intPoints(igp, :), B2, de2, vort2)

                ! get material type
                mat => mat_list(solel1%mat)				  ! element material data --> take side 1
                den1 = mat%density * his1%vol/SolidDetJ1   ! Current density
                den2 = mat%density * his2%vol/SolidDetJ2   ! Current density

                ! update stress by constitution law
                call Constitution(de1, vort1, den1, his1, mat, solidDetJ1)
                call Constitution(de2, vort2, den2, his2, mat, solidDetJ2)

                ! calculate average traction
                sm = 0.5d0 * (his1%SM + his2%SM)  ! Mean stress

                sig(1) = 0.5d0*(his1%SDS(1) + his2%SDS(1)) + sm
                sig(2) = 0.5d0*(his1%SDS(2) + his2%SDS(2)) + sm
                sig(3) = 0.5d0*(his1%SDS(3) + his2%SDS(3)) + sm
                sig(4) = 0.5d0*(his1%SDS(4) + his2%SDS(4))
                sig(5) = 0.5d0*(his1%SDS(5) + his2%SDS(5))
                sig(6) = 0.5d0*(his1%SDS(6) + his2%SDS(6))

                ! evaluate traction on the face
                traction(:) = 0.0
                traction(1) = sig(1) * faceNormal(1) + sig(4) * faceNormal(2) + sig(5) * faceNormal(3)
                traction(2) = sig(4) * faceNormal(1) + sig(2) * faceNormal(2) + sig(6) * faceNormal(3)
                traction(3) = sig(5) * faceNormal(1) + sig(6) * faceNormal(2) + sig(3) * faceNormal(3)
                
                ! add penalty contributions
                traction = traction + 1000.0 * (mat%Young / lface) * (jump_u - delta) ! beta = 4 in the DG formulation
                
                ! update interface history variables
                intfel%historyVars(igp)%jump_u = jump_u
                intfel%historyVars(igp)%p_avg = traction      
                
                if (tf > 0.5d0) traction = 0.0d0
                                
                ! evaluate and accumulate nodal forces due to element stress on the face
                call NodalForceFace(solel1%nNode, solel1%numNodes, -NshapeSolid1, traction, faceDetJ * weights_2D(igp))
                call NodalForceFace(solel2%nNode, solel2%numNodes,  NshapeSolid2, traction, faceDetJ * weights_2D(igp))
                
            end do ! end loop over integration points

    enddo ! end loop over elements

    end subroutine FEForceFaces
    
        subroutine UpdateDelta(tf)
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -    ! Update interface opening unknowns, delta                     -
    ! -                                                                   -
    ! -  Input                                                            -
    ! -    tf - value of time function at current time                    -
    ! -                                                                   -
    ! ---------------------------------------------------------------------
    use ParticleData
    use MaterialData
    use ElementData
    use Simulation
    use quadrature
    use constants	
    
    implicit none
    real(8), intent(in) :: tf
    integer :: ie, iface, igp, j, p, f, n, m, k
    integer :: numFaces, directionID
    
    real(8):: delta_c, sigma_c, delta_max, pn, pnorm			    ! density 
    integer :: solid1, solid2, face1, face2
    real(8):: xyzFace(3,4), xyzSolid1(3, 8), xyzSolid2(3, 8), v1(3,8), v2(3, 8)	! element nodal coordinates and velocities 
	real(8):: invJacobian1(3,3), invJacobian2(3,3)		! element inverse Jacobian matrix
	real(8):: de1(6), de2(6), vort1(3), vort2(3)	! incremental strain and spin
	! real(8):: pkxj(4,3)			! derivatives of shape function with respect to coordinates
    real(8):: faceDetJ, solidDetJ1,  solidDetJ2            ! determinant of the Jacobian matricies
    real(8), dimension(4, 3) :: quadratures_faces_3D ! quadrature points on the face in 3D natural coords of the element
    real(8), dimension(6, 24) :: B1, B2
    real(8), dimension(4) :: NshapeFace
    real(8), dimension(8) :: NshapeSolid1, NshapeSolid2
    real(8), dimension(3) :: jump_u, delta
    real(8), dimension(8, 3) :: dN_dXi
    real(8), dimension(3) :: faceNormal, traction
    real(8) :: lface, sm
    real(8), dimension(6) :: sig
    

    type(InterfaceElem),  POINTER :: intfel
    type(Element),  POINTER :: solel1, solel2
    type(Face), POINTER :: eface
	type(Material), POINTER :: mat
	type(History),  POINTER :: his1,his2

    ! loop over all interface elements 
    do ie = 1, nb_interface
        
        ! get intf element
        intfel => interface_list(ie)
                
        ! get solid IDS
        solid1 = intfel%solidID1
        solid2 = intfel%solidID2
        
        ! get face ID of the solids
        face1 = intfel%faceID1
        face2 = intfel%faceID2
        
        ! get the solid elements
        solel1 => element_list(solid1)
        solel2 => element_list(solid2)
        
        ! get solid elements coords and velocities
        !v1 = 0.0d0
        do j = 1, solel1%numNodes
            p = solel1%nNode(j)
            xyzSolid1(:,j) = Pos(:,p)		! Element nodal coordinates
            !v1(:,j) = Vel(:,p)		        ! Element nodal velocities
        end do
        
        !v2 = 0.0d0
        do j = 1, solel2%numNodes
            p = solel2%nNode(j)
            xyzSolid2(:,j) = Pos(:,p)		! Element nodal coordinates
            !v2(:,j) = Vel(:,p)		        ! Element nodal velocities
        end do
        
        ! loop over all faces
        do igp = 1, intfel%numIntPoints
            
            ! delta = 0.0d0
            !delta = intfel%historyVars(igp)%delta
            
            ! skip if delta_max > delta_c
            
            ! get the face nodes of the first face
        do j = 1, solel1%faces(face1)%numNodes
            p = solel1%faces(face1)%nNode(j)
                xyzFace(:,j) = Pos(:,p)		! Element nodal coordinates
        end do

        ! evluate face length
            lface = sqrt(sum((xyzFace(:, 1) - xyzFace(:, 3))**2))
        
            ! evaluate the face det(J) for surface numerical integration
            call faceJacobian(xyzFace, quadratures_2D(igp,:), solel1%faces(face1)%numNodes, faceDetJ, NshapeFace, faceNormal)
            
            !do igp = 1, eface%numIntPoints
            his1 => solel1%faces(face1)%historyVars(igp) !history_list(el%nHistory)	! element history variables
            his2 => solel2%faces(face2)%historyVars(igp)
                
            ! the element is severely distorted
            if (faceDetJ .le. 0.0d0) then
                write(*, 100) face1, ie
                100 format(1x, "*** Error *** Interface Jacobian determinant of face ", i6, " in element " , i6, " is negative ! ")
                    stop
            end if


            ! evaluate the increament strain and spin
                !call Strain_old(DTk, v, invJacobian, de, vort, pkxj)
            call solidShapeFunctions(solel1%faces(face1)%intPoints(igp, :), solel1%numNodes, NshapeSolid1, dN_dXi)
            call solidShapeFunctions(solel2%faces(face2)%intPoints(igp, :), solel2%numNodes, NshapeSolid2, dN_dXi)

            ! evaluate the jump in displacement
            jump_u = matmul(xyzSolid2, NshapeSolid2) - matmul(xyzSolid1, NshapeSolid1)

                !evaluate strain
                !call Strain(DTk, solel1%numNodes, invJacobian1, v1, solel1%faces(face1)%intPoints(igp, :), B1, de1, vort1)
                !call Strain(DTk, solel2%numNodes, invJacobian2, v2, solel2%faces(face2)%intPoints(igp, :), B2, de2, vort2)

                ! get material type
                mat => mat_list(solel1%mat)				  ! element material data --> take side 1
                sigma_c =mat%sigma_c
                delta_c = mat%delta_c
                
                !den1 = mat%density * his1%vol/SolidDetJ1   ! Current density
                !den2 = mat%density * his2%vol/SolidDetJ2   ! Current density

                ! update stress by constitution law
                !call Constitution(de1, vort1, den1, his1, mat, solidDetJ1)
                !call Constitution(de2, vort2, den2, his2, mat, solidDetJ2)
                
                ! calculate average traction
            sm = 0.5d0 * (his1%SM + his2%SM)  ! Mean stress

            sig(1) = 0.5d0*(his1%SDS(1) + his2%SDS(1)) + sm
            sig(2) = 0.5d0*(his1%SDS(2) + his2%SDS(2)) + sm
            sig(3) = 0.5d0*(his1%SDS(3) + his2%SDS(3)) + sm
            sig(4) = 0.5d0*(his1%SDS(4) + his2%SDS(4))
            sig(5) = 0.5d0*(his1%SDS(5) + his2%SDS(5))
            sig(6) = 0.5d0*(his1%SDS(6) + his2%SDS(6))

                ! evaluate traction on the face
            traction(:) = 0.0
            traction(1) = sig(1) * faceNormal(1) + sig(4) * faceNormal(2) + sig(5) * faceNormal(3)
            traction(2) = sig(4) * faceNormal(1) + sig(2) * faceNormal(2) + sig(6) * faceNormal(3)
            traction(3) = sig(5) * faceNormal(1) + sig(6) * faceNormal(2) + sig(3) * faceNormal(3)
            
            ! add penalty contributions
                traction = traction + 4.0 * (mat%Young / lface) * jump_u ! beta = 4 in the DG formulation
            
                ! update interface history variables
            intfel%historyVars(igp)%jump_u = jump_u
                intfel%historyVars(igp)%p_avg = traction     

                pn = DOT_PRODUCT(jump_u, faceNormal)
                if (pn <= 0.0d0)  traction = traction - pn * faceNormal
            
                pnorm = sqrt(sum(traction**2))
            
                if (pnorm >= sigma_c .or. intfel%historyVars(igp)%d_max > delta_c) then
                    intfel%historyVars(igp)%delta = traction / (4.0 * (mat%Young / lface))
                    intfel%historyVars(igp)%d_max = 10 ! * delta_c
                end if

        end do ! end loop over integration points

    enddo ! end loop over elements

    end subroutine UpdateDelta
    
    
    subroutine evalDeltaStar(pTraction, d_max, sigma_max, eta, delta_c, sigma_c, deltaStar)
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -    ! Evaluate the solution to opening deltaStar by solving        -
    ! -      the local non-differentiable minimization problem            -
    ! -                                                                   -
    ! -  Input:                                                           -
    ! -    pTraction - positive  part of DG traction                      -
	! -    d_max     - maximum effective opening displacement             -
    ! -    eta       - DG penalty number                                  -
    ! -    delta_c   - delta_c in the cohesive model                      -
    ! -    sigma_c   - sigma_c in the cohesive model                      -
    !    Output:                                                          -
    !-     deltaStar - optimal opening displacement                       -
    ! -                                                                   -
    ! --------------------------------------------------------------------- 
        use constants	    
        implicit none
        real(8), dimension(3), intent(in) :: pTraction
        real(8), dimension(3), intent(out) :: deltaStar
        real(8), intent(in) :: d_max, sigma_max, eta, delta_c, sigma_c
        
        
        ! local variables
        real(8) :: pNorm, pMax, Hc, dStar
        
        ! evaluate Hc
        Hc = sigma_c / delta_c
        
        ! check convexity of the objective
        if (eta <= Hc) then
            write(*, 100) eta/Hc
            100 format(1x, "*** Error *** Stability number eta is not large enough: Eta/Hc = ", f8.5, " must be > 1.0")
            stop
        end if
        
        ! evaluate norm of traction vector
        pNorm = sqrt(sum(pTraction**2))
        
        ! evaluate pMax
        pMax = sigma_max + eta * d_max
        
        if (pNorm < sigma_c .and. d_max <= EPSILON) then ! case of no fracture  (.true.) then !
            deltaStar = 0.0d0 !pTraction / eta
        else if (pNorm >= sigma_c .or. d_max > EPSILON) then
            if (pNorm <= pMax .and. d_max <= delta_c) then
                dStar = pNorm / (eta + sigma_max / d_max)
            else if (pNorm > pMax .and. pNorm <= eta*delta_c .and. d_max <= delta_c) then
                dStar = (pNorm - sigma_c) / (eta - Hc)
            else if (pNorm >= eta * delta_c .or. d_max > delta_c) then
                dStar = pNorm / eta
            else
                write(*, 101) pNorm/pMax, pNorm/(eta*delta_c)
101             format(1x, "*** Error *** Invalid traction norm: pNorm/pMax = ", f8.4, " and pNorm/eta*dc = ", f8.4)
                stop
            end if
            deltaStar = dStar * pTraction / pNorm
        else
            write(*, 102) 
102         format(1x, "*** Error *** Calculation of DeltaStar failed !")
            stop
        endif

        ! evaluate deltaStar
        

    end subroutine evalDeltaStar
    
    
 subroutine faceJacobian(nodalCoords, natCoord, numNode, detJ, Nshape, normal)
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -    ! Compute jocobian and its determinant on the element face     -
    ! -                                                                   -
    ! -  Input                                                            -
    ! -    nodalCoords - nodal coords of the face                         -
    ! -    natCoord    - natural coordinates on 2D face                   -
	! -    numNode     - number of nodes in the face                      -
    ! -    detJ        - determinant of the jacobian                      -
    ! -    normal      - normal vector to the face                        -
    ! -    Nshape      - shape functions                                  -
    ! -                                                                   -
    ! --------------------------------------------------------------------- 
  implicit none
  real(8), dimension(3, numNode), intent(in) :: nodalCoords  ! Coordinates of the 4 nodes (x, y, z)
  real(8), dimension(2), intent(in)   :: natCoord     ! Natural coordinates (xi, eta)
  
  integer, intent (in) :: numNode                     ! Number of nodes in the face
  real(8), intent(out) :: detJ                        ! Determinant of the Jacobian matrix
  real(8), intent(out), dimension(3) :: normal        ! Normal vector to the face
  real(8), dimension(numNode), intent(out) :: Nshape  ! Shape functions
  
  real(8), dimension(3) :: v1, v2, crossProd
  real(8), dimension(numNode,2) :: dN_dXi             ! Shape function derivatives w.r.t. natural coordinates
  ! real(8), dimension(numNode) :: N ! Shape functions
  real(8), dimension(3,2) :: jacobian     ! Jacobian matrix (3x2)
  integer :: i, j
  
  ! Initialize Jacobian to zero
  ! jacobian = 0.0d0
  detJ = 0.0d0
  
  call faceShapeFunctions(natCoord, numNode, Nshape, dN_dXi)
  
  ! Compute the Jacobian matrix
  do i = 1, 3
      do j = 1, 2
          jacobian(i,j) = sum(dN_dXi(:,j) * nodalCoords(i, :))
      end do
  end do
  
  ! Extract column vectors
  v1 = jacobian(:,1)
  v2 = jacobian(:,2)
  
  ! Compute cross product
  crossProd(1) = v1(2)*v2(3) - v1(3)*v2(2)
  crossProd(2) = v1(3)*v2(1) - v1(1)*v2(3)
  crossProd(3) = v1(1)*v2(2) - v1(2)*v2(1)
  
  ! Compute norm of the cross product
  detJ = sqrt(sum(crossProd**2))
  
  ! Compute the unit normal vector
  normal(1) = crossProd(1) / detJ
  normal(2) = crossProd(2) / detJ
  normal(3) = crossProd(3) / detJ
  
  ! calculate the determinant of the Jacobian matrix
  ! call calculate_determinant(jacobian, detJ)

 end subroutine faceJacobian

    
 subroutine faceShapeFunctions(natCoord, numNode, shapeFunc, dN_dXi)
   ! ---------------------------------------------------------------------
   ! -                                                                   -
   ! -  Purpose                                                          -
   ! -    ! Compute shape functions on the element face                  -
   ! -                                                                   -
   ! -  Input                                                            -
   ! -    natCoord    - natural coordinates on 2D face                   -
   ! -    numNode     - number of nodes in the face                      -
   ! -    N           - shape functions                                  -
   ! -    dN_dXi      - shape function derivatives                       -
   ! -                                                                   -
   ! ---------------------------------------------------------------------
  implicit none
  real(8), dimension(2), intent(in) :: natCoord ! Natural coordinates
  integer, intent(in) :: numNode ! Number of nodes in the face
  real(8), dimension(numNode), intent(out) :: shapeFunc ! Shape functions
  real(8), dimension(numNode, 2), intent(out) :: dN_dXi ! Shape function derivatives w.r.t. natural coordinates
  real(8):: r, s, t ! natural coordinates
  
  r = natCoord(1)
  s = natCoord(2)
  
  shapeFunc = 0.0d0
  dN_dXi = 0.0d0
  
  ! calculate the shape functions
  shapeFunc(1) = 0.25d0*(1.0d0 - r)*(1.0d0 - s)
  shapeFunc(2) = 0.25d0*(1.0d0 + r)*(1.0d0 - s)
  shapeFunc(3) = 0.25d0*(1.0d0 + r)*(1.0d0 + s)
  shapeFunc(4) = 0.25d0*(1.0d0 - r)*(1.0d0 + s)
  
  ! Calculate the shape function derivatives w.r.t. natural coordinates
  dN_dXi(1,:) = 0.25d0*(/ -(1.0d0 - s), -(1.0d0 - r) /)
  dN_dXi(2,:) = 0.25d0*(/  (1.0d0 - s), -(1.0d0 + r) /)
  dN_dXi(3,:) = 0.25d0*(/  (1.0d0 + s),  (1.0d0 + r) /)
  dN_dXi(4,:) = 0.25d0*(/ -(1.0d0 + s),  (1.0d0 - r) /)

    end subroutine faceShapeFunctions
    
    
    subroutine solidShapeFunctions(natCoord, numNode, shapeFunc, dN_dXi)
   ! ---------------------------------------------------------------------
   ! -                                                                   -
   ! -  Purpose                                                          -
   ! -    ! Compute shape functions in a solid element                   -
   ! -                                                                   -
   ! -  Input                                                            -
   ! -    natCoord    - natural coordinates on 2D face                   -
   ! -    numNode     - number of nodes in the face                      -
   ! -    N           - shape functions                                  -
   ! -    dN_dXi      - shape function derivatives                       -
   ! -                                                                   -
   ! ---------------------------------------------------------------------
  implicit none
  real(8), dimension(3), intent(in) :: natCoord ! Natural coordinates
  ! real(8), dimension(numNode, 2), intent(out) :: shapeFunctions ! Shape functions
  integer, intent(in) :: numNode ! Number of nodes in the face
  real(8), dimension(numNode, 3), intent(out) :: dN_dXi ! Shape function derivatives w.r.t. natural coordinates
  real(8), dimension(numNode), intent(out) :: shapeFunc ! Shape functions
  real(8):: r, s, t ! natural coordinates
  
  r = natCoord(1)
  s = natCoord(2)
  t = natCoord(3)
  
  shapeFunc = 0.0d0
  dN_dXi = 0.0d0
  
  ! calculate the shape functions
  shapeFunc(1) = 0.125d0*(1.0d0 - r)*(1.0d0 - s)*(1.0d0 - t)
  shapeFunc(2) = 0.125d0*(1.0d0 + r)*(1.0d0 - s)*(1.0d0 - t)
  shapeFunc(3) = 0.125d0*(1.0d0 + r)*(1.0d0 + s)*(1.0d0 - t)
  shapeFunc(4) = 0.125d0*(1.0d0 - r)*(1.0d0 + s)*(1.0d0 - t)
  shapeFunc(5) = 0.125d0*(1.0d0 - r)*(1.0d0 - s)*(1.0d0 + t)
  shapeFunc(6) = 0.125d0*(1.0d0 + r)*(1.0d0 - s)*(1.0d0 + t)
  shapeFunc(7) = 0.125d0*(1.0d0 + r)*(1.0d0 + s)*(1.0d0 + t)
  shapeFunc(8) = 0.125d0*(1.0d0 - r)*(1.0d0 + s)*(1.0d0 + t)
  
    
  dN_dXi(:,1) = 0.125d0 * (/-(1.0d0-s)*(1.0d0-t),  (1.0d0-s)*(1.0d0-t), (1.0d0+s)*(1.0d0-t), -(1.0d0+s)*(1.0d0-t), &
				-(1.0d0-s)*(1.0d0+t), (1.0d0-s)*(1.0d0+t), (1.0d0+s)*(1.0d0+t), -(1.0d0+s)*(1.0d0+t) /)
  dN_dXi(:,2) = 0.125d0 * (/-(1-r)*(1-t),  -(1+r)*(1-t), (1+r)*(1-t), (1-r)*(1-t), &
                           -(1-r)*(1+t), -(1+r)*(1+t), (1+r)*(1+t), (1-r)*(1+t) /)
  dN_dXi(:,3) = 0.125d0 * (/  -(1-r)*(1-s), -(1+r)*(1-s), -(1+r)*(1+s), -(1-r)*(1+s),  &
							(1-r)*(1-s),  (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s) /)


 end subroutine solidShapeFunctions


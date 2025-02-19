
! ---------------------------------------------------------------------
! -                                                                   -
! -  Data input procedures                                            -
! -                                                                   -
! ---------------------------------------------------------------------

module DataIn

use ElementData
use ParticleData
use MaterialData
use DataOut
use Simulation
use FFI
use constants

contains

	subroutine InputPara()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input and initialize data                                     -
! -                                                                   -
! ---------------------------------------------------------------------
	use DFLIB ! for NARGS()
	implicit none

	if(NARGS().ne.2) then
		!stop 'Usage: mpm3d InputFileName'
        FileInp = 'input.dat'
	else
		call GETARG(1,FileInp)
    end if

	call OpenInpFile()	! Open input files
    call OpenOutFile()	! Open output files

	call GetData()	! Read data from input .dat file
    
    !call GetMesh()  ! Read mesh data from .msh file
	allocate(nElem(nb_node))	! Allocate storage for stress smoothing

    end subroutine InputPara    

	subroutine GetData()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input data from input file using FFI module                   -
! -                                                                   -
! ---------------------------------------------------------------------
	use TimeFunctionData
	implicit none

	integer key, i

    !integer, parameter:: nbkw = 21
    !character(4),parameter:: kw(nbkw) = (/  &	! list of keywords
    !    'endi', 'efep', 'nbno', 'repo', 'hour', 'time', &
    !    'endt', 'outt', 'nmat', 'mate', 'node', 'jaum', &
    !    'ncom', 'load', 'ivel', 'curv', 'elem', 'nbel', &
    !    'rigi', 'dtsc', 'velo'                          &
    !    /)

	do while(.true.)
		key = keyword(kw,nbkw)
		select case(key)

		case(1)			! end input
			exit		! Terminates execution of Do loop
            !print *, 'Exit'
			
		case(2)			! epef (title)
			call GetString(Title)
            !print *, 'title'
			
		case(3)			! nbno
			nb_node = GetInt()
			print *, 'nb_node = ',nb_node

			allocate(Acc(nDim, nb_node))
			allocate(Pos(nDim, nb_node))
			allocate(Vel(nDim, nb_node))
			allocate(Fp(nDim, nb_node))
			allocate(Mp(nb_node))
            allocate(VelBCs(nb_node))

			Vel = 0.0d0
			Fp  = 0.0d0
			Mp  = 0.0d0
            !print *, 'getnode'

		case(4)			! RepTimeStep
			RepTimeStep = GetReal()
            !print *, 'time step'

		case(5)			! hourglass
			HourGlass%method = GetInt()
			HourGlass%Qhg = GetReal()
            !print *, 'hour glass'
		
		case(6)		! timefunction
			LenTMFunction = GetInt()
			allocate(TMFunction(LenTMFunction))
			do i = 1, LenTMFunction
				TMFunction(i)%time  = GetReal()
				TMFunction(i)%value = GetReal()
            end do
            !print *, 'time funciton'
			
		case(7)		! endtime
			EndTime = GetReal()
			print *, 'EndTime = ',EndTime

		case(8)		! outtimestep
			OutTimeStep = GetReal()
            !print *, 'out time step'

		case(9)		! nmat
			nb_mat = GetInt()
			print *, 'nb_mat = ',nb_mat

			allocate(mat_list(nb_mat))

		case(10)		! material
			call SetMaterial()
            !print *, 'material'

		case(11)		! node
			call SetNode()
            !print *, 'set node'

		case(12)		! Jaum
			Jaum = SetOnOff()
            !print *, 'set on-off'

		case(13)		! nb_comp
			nb_comp = GetInt()
			print *, 'nb_comp = ', nb_comp

			if (nb_node.eq.0)   &
			    stop '*** Error *** nb_node must be defined first!'

			allocate(CompLen(nb_comp))
			allocate(CompMember(nb_comp,nb_node))
			CompLen = 0
			CompMember = 0
            
		case(14)		! load
			call SetLoad()

		case(15)		! initial velocity
			call SetInitialVelocity()
            
		case(21)		! velocity BC 
			call SetVelocity()            

		case(16)		! curve
			nCurves = nCurves + 1

			if (nCurves .gt. MaxCurves)   &
				stop '*** Error *** Too many curves defined !'

			CurveOption(nCurves) = SetResOption()

			if(nb_word.gt.0) then
				CurvePoint(nCurves) = GetInt()
			else 
				CurvePoint(nCurves) = 1	  ! Default curve point
			end if

        case(17)        ! Read elements
            !print *, 'start set element'
			call ElementIn()
            !print *, 'set element'

		case(18)		! nbel
			nb_element = GetInt()
			print *, 'nb_element = ', nb_element

			allocate(element_list(nb_element))

		case(19)		! define rigid plane
			plane%nDir = GetInt()
			plane%coor = GetReal()

		case(20)		! DTSCale
			DTScale = GetReal()
        
        !case(21)		! component mesh
            !call GetMesName()
            !iComp = 
            !MeshName = 
            !GetString(CompMeshName(GetInt()))

		case default	! error
			stop 'error GetData'
			
		end select
	end do
		
    end subroutine GetData
    
!    subroutine GetMesName()
!! ---------------------------------------------------------------------
!! -                                                                   -
!! -  Purpose                                                          -
!! -     reads the name of the file for a mesh file                    -
!! -                                                                   -
!! ---------------------------------------------------------------------    
!    implicit none
!    
!    integer i
!    character(len = MAX_STRING) meshName
!
!    nb_mesh = GetInt()
!	if (nb_mesh.eq.0)   then
!            call ErrorMsg()
!            print *, '*** Error *** Number of mesh files cannot be zero.'
!            stop
!    end if
!    allocate(MeshFileName(nb_mesh))
!    MeshFileName = ''
!    
!    i = 0
!    do while(i.lt.nb_mesh)
!        i = i + 1
!        if(i.gt.nb_mesh) then
!            call ErrorMsg()
!            print *, '*** Error *** Too many mesh files listed'
!            print *, 'required : ',nb_mat
!            print *, 'inputed  : ',i
!            stop
!        end if
!
!        call GetString(meshName)
!        MeshFileName(i) = meshName
!    end do
!
!
!    end subroutine GetMesName

	subroutine SetMaterial()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input and initialize material data                            -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	integer i, t
	real(8) :: E, mu, rho

	integer,parameter:: nbkw = 8
	character(4),parameter:: kw(8) = (/'elas','pla1','pla2','john', &
                                       'elsf','plf1','plf2','jhnf' /)

	if (nb_mat.eq.0)   &
		stop '*** Error *** nb_mat must be defined in advance!'
		
	i = 0
	do while(i.lt.nb_mat)
		i = GetInt()
		if(i.gt.nb_mat) then
			call ErrorMsg()
			print *, '*** Error *** Too many material sets'
			print *, 'required : ',nb_mat
			print *, 'inputed  : ',i
			stop 
		end if

		t = KeyWord(kw,nbkw)

		select case(t)

		case(1)	! elastic
			mat_list(i)%MatType = 1
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()
            
		case(5)	! elsf: elastic with fracture
			mat_list(i)%MatType = 5
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()       
            mat_list(i)%sigma_c = GetReal()
			mat_list(i)%G_c = GetReal()   
            mat_list(i)%delta_c = 2.0 * mat_list(i)%G_c / mat_list(i)%sigma_c

		case(2) ! pla1: elastic-perfectly plastic
			mat_list(i)%MatType = 2
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()
			mat_list(i)%Yield0 = GetReal()
            
		case(6) ! plf1: elastic-perfectly plastic with fracture
			mat_list(i)%MatType = 6
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()
			mat_list(i)%Yield0 = GetReal()     
            mat_list(i)%sigma_c = GetReal()
			mat_list(i)%G_c = GetReal()          
            mat_list(i)%delta_c = 2.0 * mat_list(i)%G_c / mat_list(i)%sigma_c            

		case(3) ! pla2: isotropic hardening
			mat_list(i)%MatType = 3
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()
			mat_list(i)%Yield0 = GetReal()
			mat_list(i)%TangMod = GetReal()
            
		case(7) ! plf2: isotropic hardening with fracture
			mat_list(i)%MatType = 7
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()
			mat_list(i)%Yield0 = GetReal()
			mat_list(i)%TangMod = GetReal()  
            mat_list(i)%sigma_c = GetReal()
			mat_list(i)%G_c = GetReal()     
            mat_list(i)%delta_c = 2.0 * mat_list(i)%G_c / mat_list(i)%sigma_c

		case(4) ! john: johnson-cook
			mat_list(i)%MatType = 4
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()
			mat_list(i)%Yield0 = GetReal()
			mat_list(i)%B_jnc = GetReal()
			mat_list(i)%n_jnc = GetReal()
			mat_list(i)%C_jnc = GetReal()
            
		case(8) ! jhnf: johnson-cook with fracture
			mat_list(i)%MatType = 8
			mat_list(i)%Density = GetReal()
			mat_list(i)%Young = GetReal()
			mat_list(i)%Poisson = GetReal()
			mat_list(i)%Yield0 = GetReal()
			mat_list(i)%B_jnc = GetReal()
			mat_list(i)%n_jnc = GetReal()
			mat_list(i)%C_jnc = GetReal()      
            mat_list(i)%sigma_c = GetReal()
			mat_list(i)%G_c = GetReal()             
            mat_list(i)%delta_c = 2.0 * mat_list(i)%G_c / mat_list(i)%sigma_c

		case default
			call ErrorMsg()
			stop '*** Error *** Invalid Material Type!'

		end select

		E  = mat_list(i)%Young
		mu = mat_list(i)%Poisson
		rho= mat_list(i)%Density
	
		mat_list(i)%WaveSpeed = sqrt(E*(1-mu)/((1+mu)*(1-2*mu)*rho))

	end do

	end subroutine SetMaterial


	subroutine SetNode()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input and initialize Node data                                -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	integer i, j, p, comp, icell

	if (nb_node.eq.0)   &
		stop '*** Error *** nbmp must be defined in advance!'

	i = 0
	do while(i.ne.nb_node)
		
		i = GetInt()

		if(i.gt.nb_node) then
			call ErrorMsg()
			print *, '*** Error *** Too many nodes'
			print *, 'required : ',nb_node
			print *, 'inputed  : ',i
			stop 
		end if

		comp = GetInt()

		if(comp.gt.0 .and. comp.le.nb_comp) then
			CompLen(comp) = CompLen(comp) + 1
			CompMember(comp,CompLen(comp)) = i
		else if(comp.gt.nb_comp) then
			call ErrorMsg()
			stop '*** Error *** component number greater than nb_comp'
		end if

		do j = 1, nDim
			Pos(j,i) = GetReal()
		end do

	end do

	end subroutine SetNode


	subroutine ElementIn()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input and initialize element data                             -
! -     calculate particle mass and set cutoff mass                   -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	integer i, j, k, nbEle, nh, m, p
    integer numNodes

	real(8):: vol, mass
	real(8):: xyz(3,8), ijac(3,3)
	real(8):: xms = 0.0
    
    !integer,parameter:: nbkw = 2
	!character(4),parameter:: kw(2) = (/'hexa','tetr' /)
    
    print *, 'start set element', nb_element

	if (nb_element .eq. 0)   &
		stop '*** Error *** nbel must be defined in advance!'


	do nbEle = 1, nb_element
        
        !print *, 'element', nbEle
			
		i = GetInt()

		if(i.gt.nb_element) then
			call ErrorMsg()
			print *, '*** Error *** Too many elements'
			print *, 'required : ',nb_element
			print *, 'inputed  : ',i
			stop 
		end if

		element_list(i)%mat = GetInt()
        

		! k = KeyWord(kw,nbkw)
        k = GetInt() ! element type
        select case (k)
            case(1) ! 4-node tetrahedron
                numNodes = 4
            case(2) ! 10-node tetrahedron
                numNodes = 10
            case(3) ! 8-node hexahedral
                numNodes = 8                
		case default	! error
			call ErrorMsg()
			stop 'error: undefined element type'
		
		end select
                
        element_list(i)%numNodes = numNodes

		do j = 1, numNodes
			element_list(i)%nNode(j) = GetInt()
		enddo

    end do

!	Calculte the particle masses and volumes
    !do i = 1, nb_element
    !	do j = 1, numNodes
    !		p = element_list(i)%nNode(j)
    !		xyz(:,j) = Pos(:,p)
    !	end do
    !
    !	call SolidJacobian_old(xyz, vol, ijac)
    !
    !	history_list(element_list(i)%nHistory)%VOL = vol  ! initial element volume
    !	mass = mat_list(element_list(i)%mat)%Density * vol
    !
    !	do j = 1, element_list(i)%numNodes
    !		p = element_list(i)%nNode(j)
    !		Mp(p) = Mp(p) + mass/element_list(i)%numNodes
    !	enddo
    !
    !end do

	end subroutine ElementIn


	subroutine SetLoad()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input and initialize external load                            -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	integer k, inode, icomp, cpl, i, j
	real(8):: ff(nDim)

	integer,parameter:: nbkw = 3
	character(4),parameter:: kw(3) = (/'endl','node','comp'/)

	if (nb_node.eq.0)   &
		stop '*** Error *** nbmp must be defined !'

	do while(.true.)
		k = keyword(kw,nbkw)
		select case(k)
		case(1)	! endload
			exit

		case(2)	! by node
			inode = GetInt()
			do j = 1, nDim
				Fp(j, inode) = GetReal()
			end do

		case(3) ! by component
			icomp = GetInt()	! component number

			do j = 1, nDim
				ff(j) = GetReal()
			end do

			cpl = CompLen(icomp)! component length
			do i = 1, cpl
				inode = CompMember(icomp,i)	! node number of this component
				Fp(:,inode) = ff
			end do
				
		case default	! error
			call ErrorMsg()
			stop 'error GetLoad'
		
		end select

	end do

	end subroutine SetLoad


	subroutine SetInitialVelocity()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input and initialize velocities                               -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	integer k, inode, icomp, cpl, i, j
	real(8):: vxp(nDim)

	integer,parameter:: nbkw = 3
	character(4),parameter:: kw(3) = (/'endv','node','comp'/)

	if (nb_node.eq.0)   &
		stop '*** Error *** nbmp must be defined !'

	do while(.true.)
		k = keyword(kw,nbkw)
		select case(k)
		case(1)	! endload
			exit

		case(2)	! by node
			inode = GetInt()

			do j = 1, nDim
				Vel(j, inode) = GetReal()
			end do

		case(3) ! by component
			icomp = GetInt()	! component number
			do j = 1, nDim
				vxp(j) = GetReal()
			end do

			cpl = CompLen(icomp)  ! component length
			do i = 1, cpl
				inode = CompMember(icomp,i)	! node number of this component
				Vel(:, inode) = vxp
			end do
				
		case default	! error
			call ErrorMsg()
			stop 'error GetVelocity'
		
		end select

	end do

    end subroutine SetInitialVelocity
    
    
    subroutine SetVelocity()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Input and initialize velocity BCs                             -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	integer k, inode, iDirID, nDir, icomp, cpl, i, j
	real(8):: vxp(nDim)

	integer,parameter:: nbkw = 3
	character(4),parameter:: kw(3) = (/'endv','node','comp'/)

	if (nb_node.eq.0)   &
		stop '*** Error *** nbmp must be defined !'

	do while(.true.)
		k = keyword(kw,nbkw)
		select case(k)
		case(1)	! endload
			exit

		case(2)	! by node
			inode = GetInt()
            iDirID = GetInt()

            if (iDirID .gt. nDim .OR. iDirID .lt. 0)   then
                stop '*** Error *** Direction ID must be positive and less than nDim!'
            end if
            
            if (iDirID == 0) then  ! expect 3 velocities if dirID = 0
                do j = 1, nDim
                    VelBCs(inode)%directionIDs(j) = 1
                    VelBCs(inode)%bcVals(j) = GetReal()
                end do
            else    ! expect 1 number if 0 <= dirID <= 3 
                VelBCs(inode)%directionIDs(iDirID) = 1
                VelBCs(inode)%bcVals(iDirID) = GetReal()
            end if
    

		case(3) ! by component
			icomp = GetInt()	! component number
            iDirID = GetInt()
            
            if (iDirID .gt. nDim .OR. iDirID .lt. 0)   then
                stop '*** Error *** Direction ID must be positive and less than nDim!'
            end if
            
            if (iDirID == 0) then  ! expect 3 velocities if dirID = 0, one for each direction
                do j = 1, nDim
                    vxp(j) = GetReal()
                end do

                cpl = CompLen(icomp)  ! component length
                do i = 1, cpl
                    inode = CompMember(icomp,i)	! node number of this component
                    VelBCs(inode)%directionIDs(j) = 1
                    VelBCs(inode)%bcVals(:) = vxp
                end do
            else    ! expect 1 number if 0 <= dirID <= 3
                vxp = 0.0d0
                vxp(iDirID) = GetReal()
                cpl = CompLen(icomp)  ! component length
                do i = 1, cpl
                    inode = CompMember(icomp,i)	! node number of this component
                    VelBCs(inode)%directionIDs(iDirID) = 1
                    VelBCs(inode)%bcVals(iDirID) = vxp(iDirID)
                end do
            end if
				
		case default	! error
			call ErrorMsg()
			stop 'error GetVelocity'
		
		end select

	end do

	end subroutine SetVelocity


	logical function SetOnOff()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Read and set on/off switch                                    -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none
	integer:: k

	integer,parameter:: nbkw = 2
	character(4),parameter:: kw(2) = (/'on  ','off '/)

	k = keyword(kw,nbkw)

	select case(k)

	case(1)	! on
		SetOnOff = .true.

	case(2)	! off
		SetOnOff = .false.

	case default	! error
		call ErrorMsg()
		stop 'switch should be ON/OFF'
	
	end select

	end function SetOnOff


end module DataIn

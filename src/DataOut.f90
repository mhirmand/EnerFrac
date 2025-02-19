! ---------------------------------------------------------------------
! -                                                                   -
! -  Result output procedures                                         -
! -                                                                   -
! ---------------------------------------------------------------------

module DataOut

	use, intrinsic :: iso_fortran_env, only: int32, real32

	integer, parameter :: MaxCurves = 10	! Limited 10 curves
	integer:: nCurves = 0					! Number of curves defined
	integer:: CurveOption(MaxCurves)		! Output variable for each curve
	integer:: CurvePoint(MaxCurves)			! Particle number for curve output

	integer, parameter :: MaxAnim = 10		! Limited 10 animation variables
	integer:: nAnimate = 0					! Number of animation variables defined
	integer:: AnimateOption(MaxAnim)		! Animation variables

	real(8):: OutTimeStep = 0.0		! time interval for plot data output
	real(8):: RepTimeStep = 0.0		! time interval for status report

	integer, parameter:: nbname = 5
	character(4), parameter:: OutputName(nbname) = &
	(/'seqv', 'epef','engk', 'engs', 'ener'/)

	integer:: nRes = 0
	integer:: nPts = 0

	integer, allocatable :: nElem(:)  ! Number of elements connected with a particle

contains


	subroutine OutCurve(Time)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Output result for plotting time-history curve using TecPlot      -
! -                                                                   -
! ---------------------------------------------------------------------
	use FFI
	use ElementData, only: Element, element_list
	use MaterialData, only: History !, history_list
	use Simulation, only: Title
	implicit none

	type(Element),  POINTER :: el
	type(History),  POINTER :: his

	real(8):: i, Time

	nPts = nPts + 1

	if (nPts .eq. 1) then
		write(iow2,10) Title
		write(iow2,20) ('"', OutputName(CurveOption(i)), CurvePoint(i), '"', i=1,nCurves)
		write(iow2,*)
10		format('TITLE = "', a60, '"')
20		format('VARIABLES= "Time"', 50(A2, A4, ' at', i6, A1))
	endif

	write(iow2,"(d15.6)", advance='no') Time

	do i = 1, nCurves

		el => element_list(CurvePoint(i))
		his => element_list(CurvePoint(i))%historyVars(1) ! history_list(el%nHistory)

		call OutVariable(iow2, CurveOption(i), el%mat, his)

	end do
	write(iow2,*)

	end subroutine OutCurve


	subroutine OutAnim(Time)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Output result for post processing using TecPlot                  -
! -                                                                   -
! ---------------------------------------------------------------------
	use ParticleData
	use MaterialData
	use ElementData
	use FFI, only: iow1
	use Simulation
	implicit none

	type(Element),  POINTER :: el
	type(History),  POINTER :: his

	real(8):: Time
	integer p, nh, i, j

	nRes = nRes + 1

	if (nRes .eq. 1) then
		write(iow1,10) Title
		write(iow1,20)
		write(iow1,*)

		write(iow1, 30) Time, nb_node, nb_element
	else
		write(iow1, 40) Time, nb_node, nb_element
	endif

10	format('TITLE = "', a60, '"')
20	format('VARIABLES= "X"   "Y"   "Z"  ')
30	format('ZONE T="Time =',1p, E12.4,'"', ' F=FEPOINT N=', I5,   &
		   ' E=', I5, ' ET=BRICK C=CYAN')
40	format('ZONE T="Time =',1p, E12.4,'"', ' F=FEPOINT N=', I5,   &
		   ' E=', I5, ' ET=BRICK C=CYAN D=(FECONNECT)')

	do p = 1, nb_node
		write(iow1,"(1p, 3e12.4)") Pos(:,p)
	end do

	if (nRes .eq. 1) then
		write(iow1,"(8i7)") ((element_list(i)%nNode(j), j=1,8), i=1,nb_element)
	endif

	end subroutine OutAnim


	subroutine OutVariable(iow, opt, mat, his)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Output the value of variable name[opt]                           -
! -                                                                   -
! -  Input                                                            -
! -     iow  - File unit                                              -
! -     opt  - Variable index in name list                            -
! -     pt   - Particle                                               -
! -     his  - History varialbe                                       -
! -                                                                   -
! ---------------------------------------------------------------------
	use MaterialData, only: History
	use Simulation
	implicit none

	integer, intent(in) :: iow, opt, mat
	type(History),  intent(in) :: his

	select case(opt)

		case(1) ! seqv
			write(iow,"(d15.6)", advance='no') his%seqv
		case(2) ! epeff
			write(iow,"(d15.6)", advance='no') his%epeff
		case(3) ! engk
			write(iow,"(d15.6)", advance='no') EngKinetic
		case(4) ! engs
			write(iow,"(d15.6)", advance='no') EngStrain
		case(5) ! ener
			write(iow,"(d15.6)", advance='no') EngStrain+EngKinetic
		case default
			stop "***Error*** Invalid output variable !"

	end select

	end subroutine OutVariable


	integer function SetResOption()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Read and set result output option                             -
! -                                                                   -
! ---------------------------------------------------------------------
	use FFI
	implicit none

	SetResOption = keyword(OutputName,nbname)

	if(SetResOption.le.0 .or. SetResOption.gt.nbname) then
		call ErrorMsg()
		stop 'Result option error'
    end if

    end function SetResOption
    
 ! subroutine write_vtk_file()
 ! 	use ParticleData
	!use MaterialData
	!use ElementData
	!use FFI, only: iow1
	!use Simulation
 !   
 !   character(100) :: filename
 !   integer(int32) :: i, j, unit_num, total_connectivity
 !   character(len=20) :: str_buffer
 !   integer(int32) :: vtk_cell_type
 !
 !   filename = 'output.vtk'
 !   unit_num = 100
 !   open(unit=unit_num, file=trim(filename), form='unformatted', access='stream', status='replace', convert='LITTLE_ENDIAN')
 !
 !   ! Write VTK header
 !   write(unit_num) '# vtk DataFile Version 3.0'//char(10)
 !   write(unit_num) 'VTK file generated by Fortran'//char(10)
 !   write(unit_num) 'BINARY'//char(10)
 !   write(unit_num) 'DATASET UNSTRUCTURED_GRID'//char(10)
 !
 !   ! Write points
 !   write(str_buffer, '(I20)') nb_node
 !   write(unit_num) 'POINTS '//trim(adjustl(str_buffer))//' float'//char(10)
 !   write(unit_num) real(Pos, kind=real32)
 !
 !   ! Calculate total connectivity size
 !   total_connectivity = sum(element_list%numNodes) + nb_element
 !
 !   ! Write cells
 !   write(str_buffer, '(I20)') nb_element
 !   write(unit_num) 'CELLS '//trim(adjustl(str_buffer))//' '//trim(adjustl(str_buffer))
 !   write(str_buffer, '(I20)') total_connectivity
 !   write(unit_num) trim(adjustl(str_buffer))//char(10)
 !   do i = 1, nb_element
 !     write(unit_num) int(element_list(i)%numNodes, kind=int32), (int(element_list(i)%nNode(j)-1, kind=int32), j=1,element_list(i)%numNodes)
 !   end do
 !
 !   ! Write cell types
 !   write(str_buffer, '(I20)') nb_element
 !   write(unit_num) 'CELL_TYPES '//trim(adjustl(str_buffer))//char(10)
 !   do i = 1, nb_element
 !     vtk_cell_type = get_vtk_cell_type(element_list(i)%numNodes)
 !     write(unit_num) int(vtk_cell_type, kind=int32)
 !   end do
 !
 !   ! Write point data
 !   write(str_buffer, '(I20)') nb_node
 !   write(unit_num) 'POINT_DATA '//trim(adjustl(str_buffer))//char(10)
 !   write(unit_num) 'SCALARS nodal_data float'//char(10)
 !   write(unit_num) 'LOOKUP_TABLE default'//char(10)
 !   write(unit_num) real(Vel(1,:), kind=real32)
 !
 !   ! Write cell data
 !   write(str_buffer, '(I20)') nb_element
 !   write(unit_num) 'CELL_DATA '//trim(adjustl(str_buffer))//char(10)
 !   write(unit_num) 'SCALARS element_data float'//char(10)
 !   write(unit_num) 'LOOKUP_TABLE default'//char(10)
 !   write(unit_num) real(Vel(1,1:nb_element), kind=real32)
 !
 !   close(unit_num)
 ! end subroutine write_vtk_file
  
    subroutine write_vtk_file_ascii(nStep)
    use ParticleData
	use MaterialData
	use ElementData
	use FFI, only: FileInp
	use Simulation
    
    
    character(100) :: filename, fName
    integer :: i, j, unit_num
    integer :: nStep
    
    !! get the jobname
    fName = FileInp(1:index(FileInp,".")-1)
    
    !filename = fName // "_output_" // nStep // ".vtk"
    write(filename, '(A,A,I0,A)') trim(fName), '_output_', iStep, '.vtk'
    
    unit_num = 100*nStep + 1
    open(unit=unit_num, file=trim(filename), status='replace', action='write')

    ! Write VTK header
    write(unit_num, '(a)') '# vtk DataFile Version 3.0'
    write(unit_num, '(a)') 'VTK file generated by Fortran'
    write(unit_num, '(a)') 'ASCII'
    write(unit_num, '(a)') 'DATASET UNSTRUCTURED_GRID'

    ! Write points
    write(unit_num, '(a,i0,a)') 'POINTS ', nb_node, ' float'
    do i = 1, nb_node
      write(unit_num, '(3f12.6)') Pos(:, i)
    end do

    ! Write cells
    write(unit_num, '(a,i0,1x,i0)') 'CELLS ', nb_element, sum(element_list%numNodes) + nb_element
    do i = 1, nb_element
      write(unit_num, '(*(i0,1x))') element_list(i)%numNodes, (element_list(i)%nNode(j)-1, j=1,element_list(i)%numNodes)
    end do

    ! Write cell types
    write(unit_num, '(a,i0)') 'CELL_TYPES ', nb_element
    do i = 1, nb_element
      write(unit_num, '(i0)') get_vtk_cell_type(element_list(i)%numNodes)
    end do

    ! Write point data
    write(unit_num, '(a,i0)') 'POINT_DATA ', nb_node
    write(unit_num, '(a)') 'SCALARS nodal_data float 1'
    write(unit_num, '(a)') 'LOOKUP_TABLE default'
    do i = 1, nb_node
      write(unit_num, '(f12.6)') Vel(1, i)
    end do

    ! Write cell data
    write(unit_num, '(a,i0)') 'CELL_DATA ', nb_element
    write(unit_num, '(a)') 'SCALARS element_data float 1'
    write(unit_num, '(a)') 'LOOKUP_TABLE default'
    do i = 1, nb_element
      write(unit_num, '(f12.6)') element_list(i)%historyVars(1)%sm + element_list(i)%historyVars(1)%SDS(1)
    end do

    close(unit_num)
  end subroutine write_vtk_file_ascii

  function get_vtk_cell_type(num_nodes) result(cell_type)
  integer, intent(in) :: num_nodes
  integer :: cell_type

  select case(num_nodes)
  case(2)
      cell_type = 3  ! VTK_LINE
  case(3)
      cell_type = 5  ! VTK_TRIANGLE
  case(4)
      cell_type = 10 ! VTK_TETRA
  case(8)
      cell_type = 12 ! VTK_HEXAHEDRON
      case default
      cell_type = 7  ! VTK_POLYGON (for any other number of nodes)
  end select
  end function get_vtk_cell_type



end module DataOut
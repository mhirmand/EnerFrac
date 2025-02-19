module ParticleData
	use constants
    
	integer, parameter :: nDim = MAX_DIMENSION;		! Number of dimension
	integer :: nb_node = 0				            ! number of nodes

	real(8), allocatable :: Acc(:,:)	! particle acceleration
	real(8), allocatable :: Pos(:,:)	! particle position
	real(8), allocatable :: Vel(:,:)	! particle velocity
	real(8), allocatable :: Fp(:,:)		! external load amplitude
	real(8), allocatable :: Mp(:)		! Mass of material point
    
    
    ! Type veloBC is used to describe velocity boundary conditions
	type VelocityBC
        integer:: directionIDs(MAX_DIMENSION) ! = -1
        real(8):: bcVals(MAX_DIMENSION)  ! = 0.0d0
    end type VelocityBC 
    
    type(VelocityBC), allocatable:: VelBCs(:)		! Array of velocity boundary conditions
    
    
    !	Type RIGIDPLANE is used to describe a rigid plane (wall)
	type RigidPlane
		integer:: nDir	! The unit normal vector of the rigid plane (¡À1,¡À2,¡À3)
						!  Must point towards the impacting body
		real(8):: coor  ! x(¡À1), y(¡À2) or z(¡À3) coordinate of a point on the plane
	end type RigidPlane

	type(RigidPlane) :: plane = RigidPlane(0, 0)

end module ParticleData

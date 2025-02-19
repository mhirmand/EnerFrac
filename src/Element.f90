	
    
    module ElementData
    
    use constants
    use MaterialData
    use quadrature
    
    real(8), allocatable:: quadratures_3D(:, :), weights_3D(:)
    real(8), allocatable:: quadratures_2D(:, :), weights_2D(:)
    
    !	Type EDGE is used to describe an edge surface
	type Face
        integer:: numNodes		! number of nodes in the element
		integer:: nNode(MAX_ELE_NODES)		!  Element connectivity data
        integer:: interfaceID
        integer:: directionID
        integer:: numIntPoints
        real(8):: intPoints(MAX_INT_POINTS, 3)
        type(History):: historyVars(MAX_INT_POINTS) 
    end type Face  
    
!	Type ELEMENT is used to describe a element
	type Element
        integer:: numNodes		            ! number of nodes in the element
		integer:: nNode(MAX_ELE_NODES)		!  Element connectivity data
		integer:: mat		                ! material set number (1 ~ nb_mat)
        type(History):: historyVars(MAX_INT_POINTS) 
        type(Face):: faces(MAX_ELE_EDGES)  !
        integer:: numFaces
        integer:: numIntPoints
    end type Element  
    
    !	Type ELEMENT is used to describe a element
	type InterfaceElem
        integer:: numNodes		            ! number of nodes in the element
		integer:: nNode(MAX_ELE_NODES)		!  Element connectivity data
		integer:: mat		                ! material set number (1 ~ nb_mat)
        type(InterfaceHistory):: historyVars(MAX_INT_POINTS) 
        integer:: numIntPoints
        integer:: solidID1, faceID1
        integer:: solidID2, faceID2
    end type InterfaceElem  

	integer:: nb_element                                    ! number of elements
	type(Element), TARGET, allocatable:: element_list(:)    ! list of all element list
    
	real(8):: EleDistortion	! element distortion. 1 for no distortion
    
    integer nb_boundarySurface
    type(Face), TARGET, allocatable :: boundary_list(:)
    
    
    ! information for interface elements
    integer:: nb_interface	                                          ! number of interface elements
    type(InterfaceElem), TARGET, allocatable:: interface_list(:)      ! list of interface elements
    integer, allocatable :: nodeorigin(:,:)                           ! list of origin node of interface elements

!	Type RIGIDPLANE is used to describe a rigid plane (wall)
	type RigidPlane
		integer:: nDir	! The unit normal vector of the rigid plane (¡À1,¡À2,¡À3)
						!     Must point towards the impacting body
		real(8):: coor  ! x(¡À1), y(¡À2) or z(¡À3) coordinate of a point on the plane
	end type RigidPlane

	type(RigidPlane) :: plane = RigidPlane(0,0)

	integer:: nb_comp = 0							! number of components
	integer, allocatable:: CompLen(:)				! nb_comp
	integer, allocatable:: CompMember(:,:)          ! nb_comp * nb_node
    
    integer:: nb_mesh = 0
    character(len=MAX_STRING), allocatable:: MeshFileName(:)	! nb_comp

end module ElementData

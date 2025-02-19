    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Preprocess procedures                                            -
    ! -                                                                   -
    ! ---------------------------------------------------------------------

    module PreProcess

    use ElementData
    use ParticleData
    use MaterialData
    use Simulation
    use constants
    !use DataIn

    contains

    subroutine InitializeData()
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -     duplicate nodes to create interface elements in the mesh      -
    ! -                                                                   -
    ! ---------------------------------------------------------------------
    
    implicit none

    integer i, p, j, k, originalNode
    logical::isFracture = .false.
    real(8), allocatable:: Vel_new(:,:)
    type(velocityBC), allocatable:: VelBCs_new(:)
    
    real(8):: detJ, mass, vol, c, EleDistortion
	real(8):: xyz(3,8), ijac(3,3)
    integer :: nb_interfaceDOFs = 0, nb_stiffness = 0
        
    ! allocate and create integration point arrays. 
    ! 3D points (for solid elements)
    allocate(quadratures_3D(NumVolIntPoints**3, 3), weights_3D(NumVolIntPoints**3))
    call get_3D_integration_points(NumVolIntPoints, quadratures_3D, weights_3D)

    ! 2D domains (for interface elements).
    allocate(quadratures_2D(NumFaceIntPoints**2, 2), weights_2D(NumFaceIntPoints**2))
    call get_2D_integration_points(NumFaceIntPoints, quadratures_2D, weights_2D)

    do i = 1, nb_mat
        if(mat_list(i)%MatType.gt.4) then
            isFracture = .true.
            exit
        end if
    end do

    if (isFracture == .true.) then

        ! duplicate nodes and assign interface properties.
        if (element_list(1)%numNodes == 8) then
            call duplicateNodes_hexa()
        else if (element_list(1)%numNodes == 4) then
            call duplicateNodes_tetra()
        end if


        ! allocate now velocity and acceleration arrays
        allocate(Vel_new(nDim, nb_node))
        allocate(VelBCs_new(nb_node))
        Vel_new = 0.0d0
        do i = 1, nb_node
            originalNode = nodeorigin(i,1)
            if (Vel(1, originalNode) .ne. 0.0d0 .OR. Vel(2, originalNode) .ne. 0.0d0 .OR. Vel(3, originalNode) .ne. 0.0d0) then
                Vel_new(:,i) = Vel(:, originalNode)
            end if            
            velBCs_new(i)%directionIDs(:)  = velBCs(originalNode)%directionIDs(:)
            velBCs_new(i)%bcVals(:)  = velBCs(originalNode)%bcVals(:)
        end do


        ! replace the old particle velocity array with the new one
        if (allocated(Vel)) deallocate(Vel)
        call MOVE_ALLOC(Vel_new, Vel)
        
        !
        if (allocated(VelBCs)) deallocate(VelBCs)
        call MOVE_ALLOC(VelBCs_new, VelBCs)

        !
        if (allocated(Acc)) deallocate(Acc)
        allocate(Acc(nDim, nb_node))
        Acc = 0.0d0

        !
        if (allocated(Fp)) deallocate(Fp)
        allocate(Fp(nDim, nb_node))
        Fp = 0.0d0

        !
        if (allocated(Mp)) deallocate(Mp)
        allocate(Mp(nb_node))
        Mp = 0.0d0
    end if

    !Calculte the particle masses and volumes
    ! initiate the history variables

    EleDistortion = 1.0
    ! DTk  = DTk1
    DTk = 1.0d6
    
    do i = 1, nb_element
        element_list(i)%numIntPoints = NumVolIntPoints**3

        do j = 1, element_list(i)%numNodes
            p = element_list(i)%nNode(j)
            xyz(:,j) = Pos(:,p)
        end do
        
        mass = 0.0d0
        vol = 0.0d0
        do j = 1, element_list(i)%numIntPoints
            
            !call solidJacobian_old(xyz, vol, ijac)
            call solidJacobian(xyz, element_list(i)%numNodes, quadratures_3D(j,:), ijac, detJ)

            ! add up the element mass and volume
            mass = mass + mat_list(element_list(i)%mat)%Density * detJ * weights_3D(j)
            vol = vol + detJ * weights_3D(j)

            call InitHistory(element_list(i)%historyVars(j))
            element_list(i)%historyVars(j)%sig_y = mat_list(element_list(i)%mat)%Yield0
            element_list(i)%historyVars(j)%detJ = detJ ! / element_list(i)%numIntPoints  ! initial element volume
        enddo
        
        do j = 1, element_list(i)%numNodes
            p = element_list(i)%nNode(j)
            Mp(p) = Mp(p) + mass/element_list(i)%numNodes
        enddo

        do j = 1, element_list(i)%numFaces
            if (element_list(i)%faces(j)%directionID == 0) cycle
            do k = 1, element_list(i)%faces(j)%numIntPoints
                call solidJacobian(xyz, element_list(i)%numNodes, element_list(i)%faces(j)%intPoints(k, :), ijac, detJ)
                call InitHistory(element_list(i)%faces(j)%historyVars(k))
                element_list(i)%faces(j)%historyVars(k)%sig_y = mat_list(element_list(i)%mat)%Yield0
                element_list(i)%faces(j)%historyVars(k)%detJ = detJ ! initial element volume
            enddo
        enddo
        
        ! get wave speed
        c = mat_list(element_list(i)%mat)%WaveSpeed	
        
        ! update initial time step size
        call ElementTimeStep(xyz, vol, c, DTk, EleDistortion)

    end do



    end subroutine InitializeData
    
    
    
    subroutine duplicateNodes_tetra()
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -     duplicate nodes in a hexahedral mesh                          -
    ! -                                                                   -
    ! ---------------------------------------------------------------------

    

    end subroutine duplicateNodes_tetra


    subroutine duplicateNodes_hexa()
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -     duplicate nodes in a hexahedral mesh                          -
    ! -                                                                   -
    ! ---------------------------------------------------------------------
    
    ! Local variables
    integer :: newnnode, e, L, face, v, p, v1, v2, v3, v4, vv1, vv2
    integer :: interfacecount, boundarycount, i, ele, j
    integer :: facenodes(4,6)
    integer :: L1, L2, L3, L4, e2, face2
    integer, allocatable :: boundaryinfo(:,:)
    integer, dimension(4) :: otherkey, otherkey1, otherkey2, otherkey3, otherkey4
    integer, dimension(4) :: otherfacenodes
    real(8), dimension(4)  :: shapeFunc2D ! Shape functions
    real(8), dimension(4, 2) :: dN_dXi2D ! Shape function derivatives w.r.t. natural coordinates
    real(8):: xyzFace(3, 4), xyzSolid(3, 8), xbar(3)	! element nodal coordinates and velocities 
    real(8), dimension(3)  :: natCoords
    logical, allocatable :: paired_faces(:, :)
    
    type(Element), target, allocatable:: element_list_new(:)    ! new list of elements
    type(InterfaceElem), target, allocatable:: interface_list_temp(:)    ! new list of elements
    real(8), allocatable :: Pos_new(:,:)	                    ! new particle position array
    
    
    ! Dictionary-like structure for facelookup
    type :: dict_entry
        integer :: key(4)
        integer :: value(2)
        ! logical :: isPaired
    end type dict_entry
    type(dict_entry), allocatable :: facelookup(:)
    integer :: facelookup_count
    logical :: found_other_key, found_other_key1, found_other_key2, found_other_key3, found_other_key4
    
    ! Allocate and initialize temporary arrays
    allocate(nodeorigin(nb_element*8, MAX_DIMENSION))
    nodeorigin = 0
    
    ! Define facenodes (assuming standard hexahedral ordering: ... 
    !                   ... CCW when looked at the face from outside of the element)
    facenodes(:, 1) = (/1,4,3,2/) ! Bottom face: z = -1.0
    facenodes(:, 2) = (/5,6,7,8/) ! Top face   : z = +1.0
    facenodes(:, 3) = (/1,2,6,5/) ! Front face : x = +1.0
    facenodes(:, 4) = (/4,8,7,3/) ! Back face  : x = -1.0
    facenodes(:, 5) = (/1,5,8,4/) ! Left face  : y = -1.0
    facenodes(:, 6) = (/2,3,7,6/) ! Right face : y = +1.0
    
    newnnode = 0
    allocate(element_list_new(nb_element))
    ! element_list_new = 0

    do e = 1, nb_element
        do L = 1, 8
            newnnode = newnnode + 1
            v = element_list(e)%nNode(L)
            nodeorigin(newnnode,:) = (/v, e, L/)
            element_list_new(e)%nNode(L) = newnnode
            !element_list_new(e)%
        end do
        element_list_new(e)%numNodes = element_list(e)%numNodes
        element_list_new(e)%mat = element_list(e)%mat
        element_list_new(e)%numIntPoints = NumVolIntPoints**3; !element_list(e)%numIntPoints
        element_list_new(e)%numFaces = 6
    end do
    
    ! Initialize facelookup
    allocate(facelookup(6*nb_element))
    facelookup_count = 0
    
    do e = 1, nb_element
        do face = 1, 6
            v1 = element_list(e)%nNode(facenodes(1,face))
            v2 = element_list(e)%nNode(facenodes(2,face))
            v3 = element_list(e)%nNode(facenodes(3,face))
            v4 = element_list(e)%nNode(facenodes(4,face))
            
            ! Check if key already exists
            found_other_key = .false.
            do j = 1, facelookup_count
                if (all(facelookup(j)%key == (/v1, v2, v3, v4/))) then
                    found_other_key = .true.
                    stop
                end if
            end do
            
            if (.not. found_other_key) then
                facelookup_count = facelookup_count + 1
                facelookup(facelookup_count)%key = (/v1, v2, v3, v4/)
                facelookup(facelookup_count)%value = (/e, face/)
            end if
            
        end do
    end do
    
    ! Initialize interfacemesh and boundaryinfo
    allocate(interface_list_temp(6*nb_element))
    ! interface_list_temp = 0
    
    allocate(boundaryinfo(6*nb_element, 2))
    
    boundaryinfo = 0
    boundarycount = 0
    interfacecount = 0
    
    allocate(paired_faces(nb_element, 6))
    paired_faces = .false. 
    
    do e = 1, nb_element
        do face = 1, 6
            
            if (paired_faces(e, face) == .true.) cycle ! consider each face only once, skipping the one where vv2 > vv1 on the first edge of the face

            v1 = element_list(e)%nNode(facenodes(1,face))
            v2 = element_list(e)%nNode(facenodes(2,face))
            v3 = element_list(e)%nNode(facenodes(3,face))
            v4 = element_list(e)%nNode(facenodes(4,face))

            ! Check if the same face (v4,v3,v2,v1) exists in facelookup
            ! There could be four permutations of the same nodes in other faces
            otherkey1 = (/v4, v3, v2, v1/)
            otherkey2 = (/v3, v2, v1, v4/)
            otherkey3 = (/v2, v1, v4, v3/)
            otherkey4 = (/v1, v4, v3, v2/)

            found_other_key = .false.
            do j = 1, facelookup_count

                ! check for all possible other faces
                found_other_key1 = all(facelookup(j)%key == otherkey1)
                found_other_key2 = all(facelookup(j)%key == otherkey2)
                found_other_key3 = all(facelookup(j)%key == otherkey3)
                found_other_key4 = all(facelookup(j)%key == otherkey4)

                if (found_other_key1) then
                    otherkey = otherkey1
                    vv1 = otherkey1(3) ! appearance on the other face.
                    vv2 = otherkey1(4)
                    e2 = facelookup(j)%value(1)
                    face2 = facelookup(j)%value(2)
                    found_other_key = .true.
                    exit
                elseif(found_other_key2) then
                    otherkey = otherkey2
                    vv1 = otherkey2(2)
                    vv2 = otherkey2(3)
                    e2 = facelookup(j)%value(1)
                    face2 = facelookup(j)%value(2)
                    found_other_key = .true.
                    exit
                elseif(found_other_key3) then
                    otherkey = otherkey3
                    vv1 = otherkey3(1)
                    vv2 = otherkey3(2)
                    e2 = facelookup(j)%value(1)
                    face2 = facelookup(j)%value(2)
                    found_other_key = .true.
                    exit
                elseif(found_other_key4) then
                    otherkey = otherkey4
                    vv1 = otherkey4(4)
                    vv2 = otherkey4(1)
                    e2 = facelookup(j)%value(1)
                    face2 = facelookup(j)%value(2)
                    found_other_key = .true.
                    exit
                endif 

            enddo

            ! no other key was found. Face is a boundary face.
            if (.not. found_other_key) then 
                boundarycount = boundarycount + 1
                boundaryinfo(boundarycount,:) = (/e, face/)

                ! assign direction of the edge
                element_list_new(e)%faces(face)%directionID = 0  ! 0 dge is not on the interface (boundary edge)
                element_list_new(e)%faces(face)%numIntPoints = 0 ! no integration point
                element_list_new(e)%faces(face)%nNode(1:4) = (/v1, v2, v3, v4/) ! no integration point
                element_list_new(e)%faces(face)%intPoints = 0.0d0 ! shall be updated later
                element_list_new(e)%faces(face)%interfaceID = -1  ! no interface element exist on the face
                cycle
            end if
            
            ! update pairded face flag
            paired_faces(e, face) = .true.
            paired_faces(e2, face2) = .true.
            
            ! create interface element
            interfacecount = interfacecount + 1
            interface_list_temp(interfacecount)%numNodes = 8
            interface_list_temp(interfacecount)%solidID1 = e
            interface_list_temp(interfacecount)%solidID2 = e2
            interface_list_temp(interfacecount)%faceID1 = face
            interface_list_temp(interfacecount)%faceID2 = face2
            interface_list_temp(interfacecount)%mat = element_list(e)%mat
            interface_list_temp(interfacecount)%numIntPoints = 4            
            do j = 1, interface_list_temp(interfacecount)%numIntPoints
                call initInterfaceHistory(interface_list_temp(interfacecount)%historyVars(j))
            end do

            ! assign nodes of the interface element
            if (found_other_key1) then
                otherfacenodes(1) = facenodes(4, face2)
                otherfacenodes(2) = facenodes(1, face2)
                otherfacenodes(3) = facenodes(2, face2)
                otherfacenodes(4) = facenodes(3, face2)
            elseif (found_other_key2) then
                otherfacenodes(1) = facenodes(3, face2)
                otherfacenodes(2) = facenodes(2, face2)
                otherfacenodes(3) = facenodes(1, face2)
                otherfacenodes(4) = facenodes(4, face2)
            elseif (found_other_key3) then
                otherfacenodes(1) = facenodes(2, face2)
                otherfacenodes(2) = facenodes(1, face2)
                otherfacenodes(3) = facenodes(4, face2)
                otherfacenodes(4) = facenodes(3, face2)
            elseif (found_other_key4) then
                otherfacenodes(1) = facenodes(1, face2)
                otherfacenodes(2) = facenodes(4, face2)
                otherfacenodes(3) = facenodes(3, face2)
                otherfacenodes(4) = facenodes(2, face2)
            end if
            
            do i = 1, 4
                if (element_list(e).nNode(facenodes(i,face)) .ne. element_list(e2).nNode(otherfacenodes(i))) then
                    write(*, 100) 
                    100 format("*** Error *** Origin of the nodes in an interface element does not match !")
                    stop
                endif
            end do
            !@assert all(orgconectivity[e, thislocnums] .==
            !orgconectivity[e2, otherfacenodes])
            
            ! interface_list_temp(interfacecount)%nNode(1:4) = (/v1, v2, v3, v4/) !element_list_new(e)%nNode(facenodes(i),face))
            do i = 1, 4
                interface_list_temp(interfacecount)%nNode(i) = element_list_new(e)%nNode(facenodes(i,face))
                interface_list_temp(interfacecount)%nNode(i + 4) = element_list_new(e2)%nNode(otherfacenodes(i))
            end do

            ! assign nodes of the faces
            element_list_new(e)%faces(face)%numNodes = 4 
            element_list_new(e2)%faces(face2)%numNodes = 4 
            do i = 1, 4
                element_list_new(e)%faces(face)%nNode(i) = interface_list_temp(interfacecount)%nNode(i)
                element_list_new(e2)%faces(face2)%nNode(i) = interface_list_temp(interfacecount)%nNode(i + 4)
            end do
            
            ! assign direction of the faces (first or the second side of the interface)
            ! nodes on the  FIRST edge has a CCW ordering when looked at from outside of the element
            ! nodes on the SECOND edge has a CW ordering when looked at from outside of the element
            element_list_new(e)%faces(face)%directionID = -1  ! -1 edge is the negative side of the interface (i.e., the 1st side)
            element_list_new(e2)%faces(face2)%directionID = +1   ! +1 edge is the positive side of the interface (i.e., the 2nd side)
                            
            ! assign number of integration points of the faces
            element_list_new(e)%faces(face)%numIntPoints = 4 ! 4 integration points on the face
            element_list_new(e2)%faces(face2)%numIntPoints = 4 ! 4 integration points on the face
            
            ! assign interface element ID to the faces of both solid elements
            element_list_new(e)%faces(face)%interfaceID = interfacecount
            element_list_new(e2)%faces(face2)%interfaceID = interfacecount
            
            ! assign integration points in local coords of each solid element
            do i = 1, element_list_new(e)%faces(face)%numIntPoints
                
                do j = 1, element_list_new(e)%faces(face)%numNodes
                    p = interface_list_temp(interfacecount)%nNode(j)
                    xyzFace(:,j) = Pos(:, nodeorigin(p, 1))		! face nodal coordinates
                end do

                call faceShapeFunctions(quadratures_2D(i,:), element_list_new(e)%faces(face)%numNodes, shapeFunc2D, dN_dXi2D)
                xbar = matmul(xyzFace, shapeFunc2D)
                
                !   1. solid element 1:
                !   get the element nodal coordinates
                xyzSolid = 0.0d0
                do j = 1, element_list_new(e)%numNodes
                    p = element_list_new(e)%nNode(j)
                    xyzSolid(:,j) = Pos(:, nodeorigin(p, 1))    ! solid Element nodal coordinates
                end do
                !   find the natural coordinates of the integration point
                call natural_coord(xyzSolid, element_list_new(e)%numNodes, xbar, natCoords)
                element_list_new(e)%faces(face)%intPoints(i, 1:3) = natCoords
                
                !   2. solid element 2:
                xyzSolid = 0.0d0
                do j = 1, element_list_new(e2)%numNodes
                    p = element_list_new(e2)%nNode(j)
                    xyzSolid(:,j) = Pos(:, nodeorigin(p, 1))    ! solid Element nodal coordinates
                end do
                !   find the natural coordinates of the integration point
                call natural_coord(xyzSolid, element_list_new(e2)%numNodes, xbar, natCoords)
                element_list_new(e2)%faces(face2)%intPoints(i, 1:3) = natCoords
            end do        

        end do
    end do

    ! Resize interfacemesh
    if (allocated(interface_list)) deallocate(interface_list)
    allocate(interface_list(interfacecount))
    interface_list = interface_list_temp(1:interfacecount)
    !call MOVE_ALLOC(interface_list_temp(1:interfacecount), interface_list)
    nb_interface = interfacecount
    
    ! Create newcoords
    nb_node = 8*nb_element
    allocate(Pos_new(3, nb_node))
    do i = 1, nb_node
        Pos_new(:,i) = Pos(:, nodeorigin(i,1))
    end do
    
    ! Create boundary_list
    nb_boundarySurface = boundarycount
    allocate(boundary_list(boundarycount))
    do i = 1, boundarycount
        ele = boundaryinfo(i,1)
        face = boundaryinfo(i,2)
        do j = 1, 4  
            boundary_list(i)%nNode(j) = element_list_new(ele)%nNode(facenodes(j,face))
        end do
        
    end do
    
    ! replace the old element list with the new one
    if (allocated(element_list)) deallocate(element_list)
    call MOVE_ALLOC(element_list_new, element_list)
    
    ! replace the old particle position array with the new one
    if (allocated(Pos)) deallocate(Pos)
    call MOVE_ALLOC(Pos_new, Pos)
    
    ! clean memory
    deallocate(facelookup)
    deallocate(boundaryinfo)
    deallocate(paired_faces)
    deallocate(interface_list_temp)
    
    end subroutine duplicateNodes_hexa
    ! ---------------------------------------------------------------------
    ! -                                                                   -
    ! -  Purpose                                                          -
    ! -     duplicate nodes in a tetrahedral mesh                         -
    ! -                                                                   -
    ! ---------------------------------------------------------------------




    end module PreProcess
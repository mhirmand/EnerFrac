! module to define constants
    module constants
    implicit none
    
    
    integer, parameter :: MAX_ELE_NODES = 10
    integer, parameter :: MAX_DIMENSION = 3
    integer, parameter :: MAX_ELE_EDGES = 8
    integer, parameter :: MAX_INT_POINTS = 8
    integer, parameter :: NumVolIntPoints = 2 
    integer, parameter :: NumFaceIntPoints = 2 
    
    integer, parameter:: MAX_STRING = 256
    
    integer, parameter:: nbkw = 21
    character(4),parameter:: kw(nbkw) = (/  &	! list of keywords
        'endi', 'efep', 'nbno', 'repo', 'hour', 'time', &
        'endt', 'outt', 'nmat', 'mate', 'node', 'jaum', &
        'ncom', 'load', 'ivel', 'curv', 'elem', 'nbel', &
        'rigi', 'dtsc', 'velo'                          &
        /)
    real, parameter :: PI = 3.14159265359
    real, parameter :: betaDG = 25.0d0
    real, parameter :: EPSILON = 1e-10
    character(len=*), parameter :: VERSION = "1.0.0"
        
end module constants
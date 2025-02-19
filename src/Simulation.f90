module Simulation
	use constants 
    
	character(MAX_STRING):: Title		! problem title
	integer:: istep = 1		! Current time step
	real(8):: DTk  = 0.0	! Time step interval for variable time step integration
	real(8):: DTk1 = 0.0
	real(8):: DTk1Old = 1.0e6
	real(8):: DTScale = 0.9			! time step size factor (<= 1.0)
	real(8):: CurrentTime = 0.0		! Time for current step
	real(8):: EndTime = 0.0			! Time for end of solution

	real(8):: EngStrain  = 0.0	! strain energy
	real(8):: EngKinetic = 0.0	! kinetic energy

end module Simulation

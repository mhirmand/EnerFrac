	subroutine innve(jac,ijac,det_)
      
      implicit none

      real(8)::jac(3,3),ijac(3,3),det_

!Compute for det      
      det_=jac(1,1)*jac(2,2)*jac(3,3)&
          +jac(2,1)*jac(3,2)*jac(1,3)&
          +jac(3,1)*jac(1,2)*jac(2,3)&
          -jac(3,1)*jac(2,2)*jac(1,3)&
          -jac(2,1)*jac(1,2)*jac(3,3)&
          -jac(1,1)*jac(3,2)*jac(2,3)
!Compute for inverse      
     
    ijac(1,1)=jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2)
	ijac(2,1)=jac(2,3)*jac(3,1)-jac(2,1)*jac(3,3)
	ijac(3,1)=jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2)
	
	ijac(1,2)=jac(1,3)*jac(3,2)-jac(1,2)*jac(3,3)
	ijac(2,2)=jac(1,1)*jac(3,3)-jac(1,3)*jac(3,1)
	ijac(3,2)=jac(1,2)*jac(3,1)-jac(1,1)*jac(3,2)

	ijac(1,3)=jac(1,2)*jac(2,3)-jac(2,2)*jac(1,3)
	ijac(2,3)=jac(1,3)*jac(2,1)-jac(1,1)*jac(2,3)
	ijac(3,3)=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)

	ijac=ijac/det_

	return

    end 
    
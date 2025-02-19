module FFI
    
    use constants

	integer:: ior   = 11	! input file (mpm)
	integer:: iores = 13	! result output file (rst)
	integer:: iow1  = 14	! result animation output file (TecPlot)
	integer:: iow2  = 15	! result curves output file (TecPlot)
    integer:: iow3  = 16	! result animation output file (ParaView)

	character(100) FileInp, FileRes

!------------------------------------------------------------------
!		line_nb - current line number
!		nb_word - word left in one line, a pointer of reading position
!		nb_read - word read in one line
!------------------------------------------------------------------

	integer line_nb, nb_word, nb_read
	character(MAX_STRING) sss
	character(20) cmd_line(15)

contains
	
	subroutine OpenInpFile()
!------------------------------------------------------------------
!		purpose: Open input .dat file and initialize
!------------------------------------------------------------------
	implicit none
    !character(MAX_STRING) FileInp

	logical ext

	inquire (file=FileInp, exist=ext)
	if (.not. ext) then
		write(*,'("Input File ", a10, "does not exist !")') FileInp
		stop
	endif

	open(ior,file=FileInp,status='old')

	line_nb = 0
	nb_word = 0
	nb_read = 0

    end subroutine OpenInpFile
    
    
    subroutine OpenOutFile()
!------------------------------------------------------------------
!		purpose: Open output files and initialize
!------------------------------------------------------------------
	implicit none

	open(iores,file='mpm3d.rst',status='unknown')
	
	open(iow1, file = 'result_anim.dat', status = 'unknown')
	open(iow2, file = 'result_curv.dat', status = 'unknown')
    
    open(iow3, file = 'result_anim.vtk', status = 'unknown')

    end subroutine OpenOutFile


	subroutine GetString(mystring)
!------------------------------------------------------------------
!	purpose: Get a string
!------------------------------------------------------------------
	implicit none
	character(MAX_STRING) mystring
	mystring = sss
	nb_word = 0

	end subroutine GetString

	integer function GetInt()
!------------------------------------------------------------------
!	purpose: Get an integer number
!	inputs:
!	outputs:
!		GetInt	- the integer number
!------------------------------------------------------------------
	implicit none

	do while(nb_word.eq.0)
		nb_word = ReadLine()
		nb_read = 0
	end do

	if(isNumber(cmd_line(nb_read+1))) then
		read(cmd_line(nb_read+1),*) GetInt ! change string to digit
		nb_word = nb_word - 1
		nb_read = nb_read + 1
	else
		print *, '*** Warning *** Non Numeric Field ',cmd_line(nb_read+1),'on line ',line_nb
		nb_word = nb_word - 1
		nb_read = nb_read + 1
	end if

	end function GetInt

	real function GetReal()
!------------------------------------------------------------------
!	purpose: Get a real number
!	inputs:
!	outputs:
!		GetReal	- the real number
!------------------------------------------------------------------
	implicit none

	do while(nb_word.eq.0)
		nb_word = ReadLine()
		nb_read = 0
	end do

	if(isNumber(cmd_line(nb_read+1))) then
		read(cmd_line(nb_read+1),*) GetReal
		nb_word = nb_word - 1
		nb_read = nb_read + 1
	else
		print *, '*** Warning *** Non Numeric Field ',  &
				cmd_line(nb_read+1),'on line ',line_nb
		nb_word = nb_word - 1
		nb_read = nb_read + 1
	end if

	end function GetReal


	integer function KeyWord(kw,nbkw)
!------------------------------------------------------------------
!	purpose: read one line from ior file and identify the keyword
!	inputs:
!		kw		- keyword set
!		nbkw	- number of keyword
!	outputs:
!		Keyword	- position of inputed keyword
!------------------------------------------------------------------
	implicit none
	character(4) kw(*)
	integer i, nbkw
	
!	nb_word = 0
	KeyWord = -1

	do while(nb_word.eq.0)
		nb_word = ReadLine()
		nb_read = 0
	end do
	do i = 1, nbkw
		if(pcomp(cmd_line(nb_read+1),kw(i),4)) then
			KeyWord = i
			nb_word = nb_word - 1
			nb_read = nb_read + 1
			exit
		end if
	end do
	if(KeyWord.eq.-1) then
		print *, '*** Warning *** invalid keyword ',cmd_line(nb_read+1),'on line ',line_nb
	end if

	end function KeyWord

	integer function ReadLine()
!------------------------------------------------------------------
!	purpose: read one line from ior file and parse it
!	inputs:
!	outputs:
!		ReadLine	- number of word in this line
!		cmd_line() is filled with words in this line
!------------------------------------------------------------------
	implicit none
	integer i, j, lensss

!	initialize cmd_line
	cmd_line = ' '
!	read one line from ior
	read(ior,'(a)',err=901,end=902) sss
	line_nb = line_nb + 1

!	Strip horizontal tab character (Ctrl-I = ASCII Character 9)
	do i = 1,256
		if(ichar(sss(i:i)).eq.9) sss(i:i) = ' '
	end do
!	Strip comments
	do i = 1, 256
		if(sss(i:i).eq.'!') then
			sss(i:256) = ' '
			exit
		end if
	end do
!	Strip leading blanks
	sss = adjustl(sss)
	lensss = len_trim(sss)
	if(lensss.eq.0) then	! blank line or comment line
		ReadLine = 0
		return
	else
		ReadLine = 1
	end if
!	Remove extra blanks before separator
	i = 1
	do while(i.lt.lensss)
		do while(sss(i:i).eq.' ' .and. (sss(i+1:i+1).eq.' ' .or. sss(i+1:i+1).eq.','))
			sss(i:lensss-1) = sss(i+1:lensss)
			sss(lensss:lensss) = ' '
			lensss = lensss - 1
		end do
		i = i + 1
	end do
!	Remove extra blanks after ',' and parse sss to cmd_line
	j = 1
	i = 1
	do while(i.lt.lensss) !i = 1, lensss-2
		if(sss(i:i).eq.',' .and. sss(i+1:i+1).eq.' ') then
			sss(i+1:lensss-1) = sss(i+2:lensss)
			sss(lensss:lensss) = ' '
			lensss = lensss - 1
		end if
		if(sss(i+1:i+1).eq.',' .or. sss(i+1:i+1).eq.' ') then
			cmd_line(ReadLine) = sss(j:i)
			j = i + 2
			ReadLine = ReadLine + 1
		end if
		i = i + 1
	end do
	cmd_line(ReadLine) = sss(j:lensss)

!	see what we have read, write them to result file
	
	write(iores,"(a)") trim(sss)

	return

!	read error encountered
901	write(iores,"('error on reading, line ',i5)") line_nb
	stop 'error on reading'

!	EOF encountered	
902	print *, 'End of File!'
	stop 'encounter EOF when reading without seeing END'

	end function ReadLine
	
	logical function pcomp(a,b,n)
!c      * * F E A P * * A Finite Element Analysis Program
!c....  Copyright (c) 1984-2000: Robert L. Taylor
!c-----[--.----+----.----+----.-----------------------------------------]
!c      Purpose: Compare character strings for match
!c               Ignores upper/lower case differences.
!c      Inputs:
!c         a(*)   - Character string 1
!c         b(*)   - Character string 2
!c         n      - Number of characters to compare
!c      Outputs:
!c         pcomp  - Flag, true if a = b
!c-----[--.----+----.----+----.-----------------------------------------]
	implicit  none
	integer   n, inc, i, ia,ib
	character a*(*),b*(*)

	pcomp = .false.

!	Compute increment between an upper and lower case letter
	inc = ichar('A') - ichar('a')

!	Compare for match
	do i = 1,n
		ia = ichar(a(i:i))
		ib = ichar(b(i:i))
!	Test all permutations of characters for match
		if(ia.ne.ib .and. ia+inc.ne.ib .and. ia.ne.ib+inc ) return
	end do

	pcomp = .true.

	end function pcomp

	logical function isNumber(num)
!------------------------------------------------------------------
!	purpose: verify a numeric field
!------------------------------------------------------------------
	implicit none
	character(20) num
	character a

	a = num(1:1)
	isNumber = .false.
	if(a.eq.'-' .or. (ichar(a).gt.47 .and. ichar(a).lt.58)) isNumber = .true.

	end function isNumber

! -------------------------------------------------------------------
! -  Purpose                                                        -
! -     Display error message for input error                       -
! -------------------------------------------------------------------
	subroutine ErrorMsg()
	implicit none
	write(*,*) 'input file error, line',line_nb
	end subroutine ErrorMsg

end module FFI
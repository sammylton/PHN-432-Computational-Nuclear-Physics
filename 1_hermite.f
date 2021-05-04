	! fortran program to plot nth order hermite polynomials
	! There are some methods possible.
	! Method: Putting arrays in a common block whose dimensions are input by you is not a good method because
	! Method: Using normal arrays as arguments of the subroutine is not a good method because
	! 1) You'll have to use parameters for they are known at compile time to be able to input dimensions at the top, and they can't be passed to the subroutine directly(can't be common).
	! So, you'll have to write the parameters again in the subroutine, but the inputs should be at the top only ideally.
	! Method: Using allocatable arrays as arguments of the subroutine is the best because
	! 1) Since they don't need parameters and normal variable dimensions can be set common, they let you input the variable dimensions once at the top
	! 2) Dimensions as variables can be put in the arguments or in a common block as used by the subroutine.
	! 3) Since allocatables can't be in a common block, we use them in the argument of the subroutine.
	
	program hermite polynomials
	implicit real*8(a-h, o-z)
	character (len=*), parameter :: file2 = 'hermite.res' ! result file name
	character (len=*), parameter :: file3 = 'hermite.plt' ! gnuplot script name
	allocatable h(:,:), x(:)
        !______________________________________________________________
	! INPUT 
        !______________________________________________________________
	m = 401 ! no. of data points
	n = 6 ! order upto which hermite polynomials are wanted, n >= 0
	xi = -2.d0 ! leftmost point on x-axis
	xf = 2.d0 ! rightmost point on x-axis 
        !______________________________________________________________
	allocate (h(m,0:n), x(m))
	call linspace(x, xi, xf, m)
	call hermite(h, n, x, m)
        !______________________________________________________________
	! PRINTING DATA FILE
        !______________________________________________________________
	open (unit=2, file=file2, status='unknown')	
	do i = 1, m
	 write (2,*) x(i), (h(i,j), j = 0, n)
	enddo
	!______________________________________________________________
	! WRITING GNUPLOT SCRIPT and calling system gnuplot to plot must come after the data file is printed
	!______________________________________________________________
	open (unit=3, file=file3, status='unknown')
	write(3,*) "xi = ", xi, "; xf = ", xf, "; n = ", n
     +             , "; file2 = '", file2, "'"
	write(3,'(200A)') "set title 'Hermite Polynomials'"
	write(3,'(200A)') "set xlabel 'x'"
	write(3,'(200A)') "set ylabel 'H(x)'"
	write(3,'(200A)') "set xrange [xi:xf]"
	write(3,'(200A)') "set yrange []"
	write(3,'(200A)') "set grid"
	write(3,'(200A)') "set key inside"
	write(3,'(200A)') "set key right bottom box 1"
        write(3,'(200A)') "set terminal png #size 1920,1080 font" ! "" = " inside write statement
     +  , trim('"Arial,18"')   
	write(3,'(200A)') "# Use PDF if you want one plot"
     +  , trim(' per page (as opposed to a multiplot)')
	write(3,'(200A)') "# 11.7 and 8.3 gives A4 landscape in inches"
	write(3,'(200A)') "# set terminal pdf size 11.7,8.3 font" 
     +  , trim('"Arial,18"')
	write(3,'(200A)') "# Remove the file name extension"
	write(3,'(200A)') "z=strlen(file2)-4"
	write(3,'(200A)') "file2_base = file2[1:z]"
        write(3,'(200A)') "# Remember to change the extension to .pdf"
     +  , trim(' for the PDF terminal and exit gnuplot')
	write(3,'(200A)') "set output file2_base.'.png'"
	write(3,'(200A)') "# for %%G in (*.res) do gnuplot.exe -e"
     +  , trim(' "fn=''%%~nG''" hermite.plt')                          ! "" = "" inside a string?
        write(3,'(200A)') "plot for [i = 2:n+2] file2 u 1:i w l title"
     +  , trim(' sprintf(''H_{%i}(x)'', i-2)')                         ! '' = ' inside a string
	write(3,'(200A)') "# %i = int, %f = float, %d = double, etc."
	!______________________________________________________________
        call system('gnuplot -p hermite.plt') ! double inverted commas and semi-colons are a must for passing constants; -p makes the window stay; -e is for passing constants.
	end program hermite polynomials

	subroutine hermite(h, n, x, m)
	! hermite(hermite matrix x along row order along col, order, x array, no. of data points)
	! To create a hermite polynomial matrix, where x changes along the rows, and the order along columns
	! It's better to write a subroutine than a function for them because they are anyway running from 0 to nth order
	implicit real*8(a-h, o-z)
	dimension h(m,0:n), x(m)
	do i = 1, m
	 h(i,0) = 1.d0 
	 if (n .gt. 0) then
       	  h(i,1) = (2.d0)*x(i)
	  do j = 2, n
	   h(i,j) = (2.d0)*(x(i)*h(i,j-1)-(j-1)*h(i,j-2))
	  enddo
	 endif
	enddo
	return 
	end subroutine hermite
	
	subroutine linspace(x, xi, xf, n)
	! linspace(array to be created, lower limit, upper limit, no. of data points)
	implicit real*8(a-h, o-z)
	dimension x(n)
	if (n == 0) return
	if (n == 1) then
	 x(1) = xi ! initializing a point on x-axis
	 return	
	endif
	do i = 1, n 
	 x(i) = xi + (xf-xi)*(i-1)/(n-1)
	! This method ensures that the array ends at xf 'very precisely'.
	enddo
	return
	end subroutine linspace
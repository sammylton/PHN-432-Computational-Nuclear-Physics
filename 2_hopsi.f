	! fortran program to plot nth order Harmonic oscillator wave function in natural units.

	program harmonic oscillator wave function
	implicit real*8 (a-h,o-z)
	real start, finish ! to calculate time
	allocatable h(:,:), psi(:), factv(:), x(:)
	character (len=*), parameter :: file2 = 'hopsi.res' ! result file name
	character (len=*), parameter :: file3 = 'hopsi.plt' ! gnuplot script name, IMP: also present at calling system gnuplot
        !______________________________________________________________
	! INPUT 
        !______________________________________________________________
	m = 4001 ! no. of data points
	n = 170 ! the largest order required of the wave function, one can input an array of orders.
	! since (170!) is nearing the maximum limit of double precision in Fortran
	xi = -20.d0 ! leftmost point on x-axis
	xf = 20.d0 ! rightmost point on x-axis
        !______________________________________________________________
	allocate (h(m,0:n), psi(m), factv(0:n), x(m))
	call linspace(x, xi, xf, m)
	call hermite(h, n, x, m)
	call factorial(factv, n)
        !______________________________________________________________
	! PRINTING DATA FILE
        !______________________________________________________________
	open (unit=2, file=file2, status='unknown')
	call howf(psi, h, factv, x, n, m) ! needs to be called for every new order value
	do i = 1, m
	 write (2,*) x(i), psi(i)
	enddo
	!______________________________________________________________
	! WRITING GNUPLOT SCRIPT 
	!______________________________________________________________
	open (unit=3, file=file3, status='unknown')
	write(3,*) "xi = ", xi, "; xf = ", xf, "; n = ", n
     +             , "; file2 = '", file2, "'"
	write(3,'(200A)') "set title 'Harmonic Oscillator Wave Functn(s)'"
	write(3,'(200A)') "set xlabel 'x'"
	write(3,'(200A)') "set ylabel '\psi(x)'"
	write(3,'(200A)') "set xrange [xi:xf]"
	write(3,'(200A)') "set yrange []"
	write(3,'(200A)') "set grid"
	write(3,'(200A)') "set key inside"
	write(3,'(200A)') "set key right bottom box 1"
        write(3,'(200A)') "set terminal png size 1920,1080 font" ! "" = " inside write statement
     +  , trim('"Arial,18"')   
	write(3,'(200A)') "# set terminal pdf size 11.7,8.3 font" 
     +  , trim('"Arial,18"')
	write(3,'(200A)') "# Remove the file name extension"
	write(3,'(200A)') "z=strlen(file2)-4"
	write(3,'(200A)') "file2_base = file2[1:z]"
	write(3,'(200A)') "set output file2_base.'.png'"                     ! "" = "" inside a string?
        write(3,'(200A)') "plot file2 u 1:2 w l title"
     +  , trim(' sprintf(''\psi_{%i}(x)'', n)')                       
	!______________________________________________________________
	call system('gnuplot -p hopsi.plt')
	end program harmonic oscillator wave function

	subroutine factorial(factv, n)
	! factorial(array to store fact(0:n), no. upto which factorials are wanted)
	! It's better to write a subroutine for factorial because it anyway calculates the factorial of smaller nos; when factorials of many nos are needed.
	implicit real*8(a-h,o-z)
	dimension factv(0:n)
	factv(0) = 1.d0 ! note that 0! = 1 
	do i = 1, n
	 factv(i)=factv(i-1)*i
	enddo
	return
	end subroutine factorial

	subroutine hermite(h, n, x, m)
	! hermite(hermite matrix x along row order along col, order, x array, no. of data points)
	! To create a hermite polynomial matrix, where x changes along the rows, and the order along columns
	! It's better to write a subroutine than a function for them because they are anyway running from 0 to nth order; when many orders are needed.
	! subroutine linspace should be called beforehand.
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

	subroutine howf(psi, h, factv, x, n, m)
	! To create the harmonic oscillator wave functions as a function of order and x
	! howf(wave function as a function of x and order, hermite polynomial matrix, factorial array, linspace x array, order of hermite polynomial, no. of data points)
	! subroutine linspace, factorial, hermite should be called beforehand.
	implicit real*8(a-h,o-z)
	parameter (pi = 4.d0*datan(1.d0)) 
	dimension h(m,0:n), psi(m), factv(0:n), x(m)
	do i = 1, m
	 psi(i) = dexp(-0.5d0*x(i)**2.d0)*h(i,n)
     +            /dsqrt(dsqrt(pi)*2.d0**n*factv(n))
	enddo
	end subroutine howf
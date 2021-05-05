	! To establish the orthonormalization condition of the Harmonic Oscillator Wave Functions 
	! and find the expectation value of Harmonic oscillator potential and kinetic terms in natural units.	

	program orthonormalization and expectation 
	implicit real*8(a-h, o-z)
	real start, finish ! to calculate time
	allocatable h(:,:), psin1(:), psin2(:), psin2_(:), factv(:)
     +              , fn(:), fv(:), ft(:), x(:)
	!______________________________________________________________
	! INPUT
	!______________________________________________________________
	m=4001 ! the no. of data points, must be odd for Simpson's
	n1 = 170
	n2 = 170! the orders of wave functions 
	xi = -20.d0 ! leftmost point on x-axis compatible with the greatest order possible
	xf = 20.d0 ! rightmost point on x-axis compatible with the greatest order possible
	!______________________________________________________________
	n = max(n1, n2)
	allocate (h(m,0:n), psin1(m), psin2(m), psin2_(m), factv(0:n)
     +           , fn(m), fv(m), ft(m), x(m))
	call linspace(x, xi, xf, m)
	call hermite(h, n, x, m)
	call factorial(factv, n)
	call howf(psin1, h, factv, x, n1, m)
	call howf(psin2, h, factv, x, n2, m)
	k = 0 ! Integrand of the Expectation of Kinetic Energy formula is different when n2 = 0
	if (n2 .ne. 0) then
	 k = 1
	 call howf(psin2_, h, factv, x, n2-1, m)
	endif 
	!______________________________________________________________
	! INTEGRANDS OF THE ORTHONORMALIZATION AND EXPECTATION FORMULAE
	!______________________________________________________________
	do i = 1, m
	 fn(i) = psin1(i)*psin2(i)
	 fv(i) = 0.5d0*x(i)**2.d0*psin1(i)*psin2(i)
	 ft(i) = 0.5d0*psin1(i)*(k*dsqrt(8.d0*n2)*x(i)*psin2_(i)
     +           +(1.d0-x(i)**2.d0)*psin2(i))
	enddo
	!______________________________________________________________
	! INTEGRATIONS
	!______________________________________________________________
	sn = simpson(fn, xi, xf, m) ! <n1|n2>
	sv = simpson(fv, xi, xf, m) ! <n1|V|n2>
	st = simpson(ft, xi, xf, m) ! <n1|T|n2>
	!______________________________________________________________
	! PRINTING OUTPUT
	!______________________________________________________________
	print*, 'normalization =', sn, 'Expectation value of V=', sv
     +        , 'Expectation value of T=', st
	end program orthonormalization and expectation

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

	real*8 function simpson(f, xi, xf, n)
	! simpson(fn being integrated, lower lim, upper lim, no. of data points)
	implicit real*8(a-h, o-z)
	dimension f(n)
	ns = n-1 ! no. of segments of the domain for Simpson's, must be even
	se = 0.d0 ! initializing 
	so = 0.d0 ! initializing 
	do i = 2, ns, 2
	 se = se + f(i)
	 so = so + f(i+1)
	enddo
        simpson=((xf-xi)*(f(1)+(4.d0)*se+(2.d0)*so-f(n)))/(ns*(3.d0))       
	! sum of alternate terms after the first gets 4X
	! Don't talk odd/even. Talk first/last and relative to them.
	return 
	end function simpson
	
	
	
	

	
	
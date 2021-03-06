! This file is provided as part of the "projet long" for the Algebre Lineaire Numerique course
! at ENSEEIHT
! Date: 04/2015
! Authors: P. Amestoy, P. Berger, A. Buttari, Y. Diouane, S. Gratton, F.H. Rouet, E. Simon
!
! This file contains the implementation of
! - one routine for computing the eigenvalues of a symmetric matrix
! - one routine for computing the singular values of an unsymmetric matrix
!
! both are based on the subspace iteration method
!
! module m_subspace_iter
  ! implicit none
  ! private
  ! public :: subspace_iter_ev, subspace_iter_sv
! contains




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subspace iteration with Rayleigh-Ritz.
! In this case the convergence can be checked for each eigenvector separately
! and the method can be stopped when the convereged eigenvectors capture
! as much as the trace of A as requested by the "percentage" argument
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine subspace_iter_ev(a, v, d, n, l, p, percentage, maxit, eps,  &
       res_ev, it_ev, it, n_ev, ierr)
    implicit none
    !! the matrix and subspace dimensions
    integer,          intent(in)                            :: n, l
    !! the number of products per iteration
    integer,          intent(in)                            :: p
    !! the traget matrix
    double precision, dimension(n,n), intent(in)            :: a
    !! maximum # of iteration
    integer,          intent(in)                            :: maxit
    !! the number of the dominant eigenvectors to compute
    double precision, intent(in)                            :: percentage
    !! the tolerance for the stopping criterion
    double precision, intent(in)                            :: eps
    !! the starting subspace. The computed eigenvectors will be
    !! returned in this array
    double precision, dimension(n,l), intent(inout)         :: v
    !! the residuals for each each eigenvector
    double precision, dimension(l),   intent(out)           :: res_ev
    !! the n_ev dominant eigenvalues
    double precision, dimension(l),   intent(out)           :: d
    !! the number of iteration to converge for each eigenvector
    integer,          dimension(l),   intent(out)           :: it_ev
    !! the global number iteration to converge
    integer,          intent(out)                           :: it
    !! the number of converged eigenvectors
    integer,          intent(out)                           :: n_ev
    !! a flag for signaling errors
    integer,          intent(out)                           :: ierr

    !! local variables
    integer                                                 :: i, j, conv, lwork
    double precision, allocatable, dimension(:,:)           :: y
    double precision, allocatable, dimension(:,:)           :: h
    double precision, allocatable, dimension(:)             :: w
    double precision                                        :: lambda, res, tmp,  trace, eig_sum
    integer, parameter                                      :: ione=1

    ! external functions
    double precision, external                              :: dlange, dnrm2

#if defined(mex)
    character(len=80) :: string
    integer*4, external :: mexPrintf
    integer*4 :: k
#endif


    trace = 0.d0
    do i=1, n
       trace = trace+a(i,i)
    end do
    eig_sum = 0.d0
        
    ierr = 0
    !! initialization process ...
    if((percentage.gt.1d0)  .or. (percentage.lt.0d0)) then
       ierr=1
       return
    end if
    if(l.gt.n) then
       ierr=1
       return
    end if

    allocate(y(n,l), h(l,l), w(l), stat=ierr)
    if(ierr .ne. 0) return

    lwork = n*l
    
    ! a natural choice of lambda is an estimate of some norm of a.
    ! use y as a workspace
    lambda=dlange('f', n, n, a, n, y)


    ! initialize the subspace
    call random_number(v)
    call orth_basis_ip(v, n, n, l)
    
    n_ev      = 0
    res_ev    = 0d0
    d         = 0d0
    it_ev     = 0
    
    it=0
    do while( (n_ev.lt.l) .and. (it.lt.maxit) )

       it = it+1
       !!compute  y=a*v
       do i=1, p
          call dgemm('n', 'n', n, l, n, 1.d0, a, n, v, n, 0.d0, y, n)
          v=y
       end do

       !! orthogonalization using Gramm Schmidt procedure
       !! compute  v=orth(y)
       call orth_basis_ip(v, n, n, l)

       !!compute  y=a*v
       call dgemm('n','n', n, l, n, 1.d0, a, n, v, n, 0.d0, y, n)
       !!compute  h=v'*y
       call dgemm('t', 'n', l, l, n, 1.d0, v, n, y, n, 0.d0, h, l)
       ! Compute eig-decomposition of h. Use y as a workspace
       call dsyev('v', 'u', l, h, l, w, y, lwork, ierr)
       if( ierr .ne.0 )then
          goto 9999
       end if

       !! Sort in decreasing order (in-place)
       !!(we suppose that all the eigen values are positive)
       do i=1,l/2
          tmp        = w(i)
          w(i)       = w(l-i +1)
          w(l-i+1)   = tmp
          ! use the first column of y as a workspace
          y(1:l,1)   = h(:,i)
          h(:,i)     = h(:,l-i +1)
          h(:,l-i+1) = y(1:l,1)
       end do

       !! v=v*h
       y = v
       call dgemm('n', 'n', n, l, l, 1d0, y, n, h, l, 0d0, v, n)

       conv = 0
       i    = n_ev +1
       !! the larger eigenvalue will converge more swiftly than
       !! those corresponding to the smaller eigenvalue.
       !! for this reason, we test the convergence in the order
       !! i=1,2,.. and stop with the first one to fail the test
       do
          if( i .gt. l) exit
          !! compute res=norm(a*v(:,i) - v(:,i)*t(i,i),2)/lambda;
          !! use the first column of y as a workspace
          y(:,1) = v(:,i)
          call dgemv('n', n, n, 1d0, a, n, v(1,i), 1, -w(i), y, 1)
          res = dnrm2(n, y, ione)/lambda

          if(res.gt.eps) exit
          conv         = conv+1
          d(i)         = w(i)
          res_ev(i)    = res
          it_ev(i)     = it
          eig_sum      = eig_sum+w(i)
          i            = i+1
       end do

       n_ev = n_ev + conv

#if defined (mex)
       write(string,'(" IT:",i5," -- Found ",i4," eigenvalues (",f8.5," percent achieved)")')it,n_ev,eig_sum/trace
       k = mexprintf(string//achar(10))
#else
       write(*,'(" IT:",i5," -- Found ",i4," eigenvalues (",f8.5,"% achieved)",a)',advance='no')it,n_ev,eig_sum/trace,char(13)
       write(*,'(" ")')
#endif
       if( eig_sum/trace .gt. percentage) exit
       if( n_ev .ge. l) exit

    end do


    if(it .ge. maxit)then
       ierr=4
       goto 9999
    end if

9999 continue
    deallocate(y, h, w)

    return

  end subroutine subspace_iter_ev



  subroutine subspace_iter_sv_left(a, u, s, m, n, l, p, percentage, maxit, eps,   &
       res_sv, it_sv, it, n_sv, ierr)
    implicit none
    !! the subspace dimensions
    integer,          intent(in)                            :: m, n, l
    !! the number of products per iteration
    integer,          intent(in)                            :: p
    !! the traget matrix
    double precision, dimension(m,n), intent(in)            :: a
    !! maximum # of iteration
    integer,          intent(in)                            :: maxit
    !! the number of the dominant eigenvectors to compute
    double precision, intent(in)                            :: percentage
    !! the tolerance for the stopping criterion
    double precision, intent(in)                            :: eps
    !! the starting subspace. The computed eigenvectors will be
    !! returned in this array
    double precision, dimension(m,l), intent(inout)         :: u
    !! the residuals for each each eigenvector
    double precision, dimension(l),   intent(out)           :: res_sv
    !! the n_sv dominant eigenvalues
    double precision, dimension(l),   intent(out)           :: s
    !! the number of iteration to converge for each eigenvector
    integer,          dimension(l),   intent(out)           :: it_sv
    !! the global number iteration to converge
    integer,          intent(out)                           :: it
    !! the number of converged eigenvectors
    integer,          intent(out)                           :: n_sv
    !! a flag for signaling errors
    integer,          intent(out)                           :: ierr

    !! local variables
    integer                                                 :: i, j, conv, lwork
    double precision, allocatable, dimension(:,:)           :: y, z
    double precision, allocatable, dimension(:,:)           :: h
    double precision, allocatable, dimension(:)             :: w
    double precision                                        :: lambda, res, tmp, trace, eig_sum
    integer, parameter                                      :: ione=1

    ! external functions
    double precision, external                              :: dlange, dnrm2

#if defined(mex)
    character(len=80) :: string
    integer*4, external :: mexPrintf
    integer*4 :: k
#endif

    ierr = 0
    !! initialization process ...
    if((percentage.gt.1d0)  .or. (percentage.lt.0d0)) then
       ierr=1
       return
    end if
    if(l.gt.n) then
       ierr=1
       return
    end if

    allocate(y(m,l), h(l,l), w(l), z(n,l), stat=ierr)
    if(ierr .ne. 0) return

    trace = dnrm2(m*n, a(1,1), 1)
    trace = trace*trace
    eig_sum = 0.d0
    
    lwork = m*l
    
    ! a natural choice of lambda is an estimate of some norm of a.
    ! use y as a workspace
    lambda=dlange('f', m, n, a, m, y)
    lambda=lambda*lambda

    call random_number(u)
    call orth_basis_ip(u, m, m, l)

    n_sv      = 0
    res_sv    = 0d0
    s         = 0d0
    it_sv     = 0
    
    it=0
    do while( (n_sv.lt.l) .and. (it.lt.maxit) )

       it = it+1
       !!compute  y=(a*a')^p * u
       do i=1, p
          ! z=a'*u
          call dgemm('t', 'n', n, l, m, 1.d0, a, m, u, m, 0.d0, z, n)
          ! u = a*z
          call dgemm('n', 'n', m, l, n, 1.d0, a, m, z, n, 0.d0, u, m)
       end do

       !! orthogonalization using Gramm Schmidt procedure
       !! compute  v=orth(y)
       call orth_basis_ip(u, m, m, l)

       ! z=a'*u
       call dgemm('t', 'n', n, l, m, 1.d0, a, m, u, m, 0.d0, z, n)
       !!compute  y=a*z
       call dgemm('n', 'n', m, l, n, 1.d0, a, m, z, n, 0.d0, y, m)
       !!compute  h=u'*y
       call dgemm('t', 'n', l, l, m, 1.d0, u, m, y, m, 0.d0, h, l)
       ! Compute eig-decomposition of h. Use y as a workspace
       call dsyev('v', 'u', l, h, l, w, y, lwork, ierr)
       if( ierr .ne.0 )then
          goto 9999
       end if

       !! Sort in decreasing order (in-place)
       !!(we suppose that all the eigen values are positive)
       do i=1,l/2
          tmp        = w(i)
          w(i)       = w(l-i +1)
          w(l-i+1)   = tmp
          ! use the first column of y as a workspace
          y(1:l,1)   = h(:,i)
          h(:,i)     = h(:,l-i +1)
          h(:,l-i+1) = y(1:l,1)
       end do

       !! u=u*h
       y = u
       call dgemm('n', 'n', m, l, l, 1d0, y, m, h, l, 0d0, u, m)

       conv = 0
       i    = n_sv +1
       !! the larger eigenvalue will converge more swiftly than
       !! those corresponding to the smaller eigenvalue.
       !! for this reason, we test the convergence in the order
       !! i=1,2,.. and stop with the first one to fail the test
       do
          if( i .gt. l) exit
          !!compute res=norm(a*v(:,i) - v(:,i)*t(i,i),2)/lambda;
          !!--compute aux_res=a*v(:,i) - v(:,i)*t(i,i)
          !!--compute res=||aux_res||/||a||
          ! use the first column of z as a workspace
          call dgemv('t', m, n, 1.d0, a, m, u(1,i), 1, 0.d0, z, 1)
          y(:,1) = u(:,i)
          call dgemv('n', m, n, 1.d0, a, m, z, 1, -w(i), y, 1)
          res = dnrm2(m, y, ione)/lambda

          if(res.gt.eps) exit
          conv         = conv+1
          s(i)         = sqrt(w(i))
          res_sv(i)    = res
          it_sv(i)     = it
          eig_sum      = eig_sum+w(i)
          i            = i+1
       end do

       n_sv = n_sv + conv

#if defined (mex)
       write(string,'(" IT:",i5," -- Found ",i4,"  singular values (",f8.5,"percent achieved)")')it,n_sv,eig_sum/trace
       k = mexprintf(string//achar(10))
#else
       write(*,'(" IT:",i5," -- Found ",i4," singular values (",f8.5,"% achieved)",a)',advance='no')it,n_sv,eig_sum/trace,char(13)
       write(*,'(" ")')
#endif

       if( eig_sum/trace .gt. percentage) exit
       if( n_sv .ge. l) exit

    end do

    if(it .ge. maxit)then
       ierr=4
       goto 9999
    end if

9999 continue
    deallocate(y, h, w, z)

    return

  end subroutine subspace_iter_sv_left

  

  subroutine subspace_iter_sv_both(a, u, v, s, m, n, l, p, percentage, maxit, eps,   &
       res_sv, it_sv, it, n_sv, ierr)
    implicit none
    !! the subspace dimensions
    integer,          intent(in)                            :: m, n, l
    !! the number of products per iteration
    integer,          intent(in)                            :: p
    !! the traget matrix
    double precision, dimension(m,n), intent(in)            :: a
    !! maximum # of iteration
    integer,          intent(in)                            :: maxit
    !! the number of the dominant eigenvectors to compute
    double precision, intent(in)                            :: percentage
    !! the tolerance for the stopping criterion
    double precision, intent(in)                            :: eps
    !! the starting subspace for the left singular vectors. The
    !! computed left singular vectors will be returned in this array
    double precision, dimension(m,l), intent(inout)         :: u
    !! the starting subspace for the right singular vectors. The
    !! computed right singular vectors will be returned in this array
    double precision, dimension(n,l), intent(inout)         :: v
    !! the residuals for each each eigenvector
    double precision, dimension(l),   intent(out)           :: res_sv
    !! the n_sv dominant eigenvalues
    double precision, dimension(l),   intent(out)           :: s
    !! the number of iteration to converge for each eigenvector
    integer,          dimension(l),   intent(out)           :: it_sv
    !! the global number iteration to converge
    integer,          intent(out)                           :: it
    !! the number of converged eigenvectors
    integer,          intent(out)                           :: n_sv
    !! a flag for signaling errors
    integer,          intent(out)                           :: ierr

    !! local variables
    integer                                                 :: i, j
    integer, parameter                                      :: ione=1, version=2
    real(kind(1.d0)), parameter                             :: done=1.d0, dzero=0.d0
    real(kind(1.d0)), allocatable                           :: b(:,:)

    select case (version)
    case (1)
       call subspace_iter_sv_left(a, u, s, m, n, l, p, percentage, maxit, eps,   &
            res_sv, it_sv, it, n_sv, ierr)
       
       call dgemm('t', 'n', n, n_sv, m, done, a, m, u, m, dzero, v, n)
       
       do j = 1, n_sv
          v(1:n,j) = v(1:n,j)/s(j)
       end do

    case (2)
       allocate(b(n, m))
       b = transpose(a)
       
       call subspace_iter_sv_left(b, v, s, n, m, l, p, percentage, maxit, eps,   &
            res_sv, it_sv, it, n_sv, ierr)

       call dgemm('n', 'n', m, n_sv, n, done, a, m, v, n, dzero, u, m)
       
       do j = 1, n_sv
          u(1:m,j) = u(1:m,j)/s(j)
       end do
       deallocate(b)

    case (3)
       allocate(b(n, n))
       call dgemm('t', 'n', n, n, m, done, a, m, a, m, dzero, b, n)

       call subspace_iter_ev(b, v, s, n, l, p, percentage, maxit, eps,   &
            res_sv, it_sv, it, n_sv, ierr)
       
       call dgemm('n', 'n', m, n_sv, n, done, a, m, v, n, dzero, u, m)
       
       do j = 1, n_sv
          s(j) = sqrt(s(j))
          u(1:m,j) = u(1:m,j)/s(j)
       end do
       deallocate(b)
    end select
    return

  end subroutine subspace_iter_sv_both




  subroutine subspace_iter_sv_both_alternate(a, u, v, s, m, n, l, p, percentage, maxit, eps,   &
       res_sv, it_sv, it, n_sv, ierr)
    implicit none
    !! the subspace dimensions
    integer,          intent(in)                            :: m, n, l
    !! the number of products per iteration
    integer,          intent(in)                            :: p
    !! the traget matrix
    double precision, dimension(m,n), intent(in)            :: a
    !! maximum # of iteration
    integer,          intent(in)                            :: maxit
    !! the number of the dominant eigenvectors to compute
    double precision, intent(in)                            :: percentage
    !! the tolerance for the stopping criterion
    double precision, intent(in)                            :: eps
    !! the starting subspace for the left singular vectors. The
    !! computed left singular vectors will be returned in this array
    double precision, dimension(m,l), intent(inout)         :: u
    !! the starting subspace for the right singular vectors. The
    !! computed right singular vectors will be returned in this array
    double precision, dimension(n,l), intent(inout)         :: v
    !! the residuals for each each eigenvector
    double precision, dimension(l),   intent(out)           :: res_sv
    !! the n_sv dominant eigenvalues
    double precision, dimension(l),   intent(out)           :: s
    !! the number of iteration to converge for each eigenvector
    integer,          dimension(l),   intent(out)           :: it_sv
    !! the global number iteration to converge
    integer,          intent(out)                           :: it
    !! the number of converged eigenvectors
    integer,          intent(out)                           :: n_sv
    !! a flag for signaling errors
    integer,          intent(out)                           :: ierr


    ! external functions
    double precision, external                              :: dlange, dnrm2
	
	!Variables locales
	integer													:: currentS, i, j, lwork, conv
	double precision										:: traceAM, traceAN, normeA, res, tmp, percentReached
	double precision, allocatable, dimension(:,:)           :: y,uout, matM,matN
	double precision, allocatable, dimension(:,:)           :: z,vout
    double precision, allocatable, dimension(:,:)           :: h
    double precision, allocatable, dimension(:)             :: w
    
    ! initialization of n_sv - mandatory for the call with mex
    n_sv=0
    
    ierr = 0
    !! initialization process ...
    if((percentage.gt.1d0)  .or. (percentage.lt.0d0)) then
       ierr=1
       return
    end if
    if(l.gt.n) then
       ierr=1
       return
    end if

   !Initialisation
	it = 0
	currentS = 1
        traceAM = 0
        traceAN = 0
	n_sv = 0
	percentReached = 0
	!Autres intialisations 
    
    !Allocation de mémoire!
	allocate(uout(m,l),vout(n,l),y(m,l),z(n,l), h(l,l), w(l),matM(m,m), matN(n,n), stat=ierr)
    if(ierr .ne. 0) return
    
    ! initialize the subspace

    call random_number(vout)
    call random_number(uout)
        
    call orth_basis_ip(vout, n,n, l)
    call orth_basis_ip(uout, m,m, l)
    
    lwork = m*l

     !Calcul de la trace de A!
          call dgemm('n', 't', m, m, n, 1.d0, a, m, a, m, 0.d0, matM, m)
          traceAM = dlange('f',m,m,MatM,m,y)
          call dgemm('t', 'n', n, n, m, 1.d0, a, m, a, m, 0.d0, matN, n)
          traceAN = dlange('f',n,n,MatN,n,z)

	! utilisation de u pour stockage!
	do while( n_sv.lt.l .and. it <= maxit )
		if (currentS == 1) then
			do i=1,p
				if ( mod(i,2) == 0) then
					call dgemm('n', 'n', m, l, n, 1.d0, a, m, vout, n, 0.d0, uout, m)
				else
				        call dgemm('t', 'n', n, l, m, 1.d0, a, m, uout, m, 0.d0, vout, n)
				end if
			end do
		
		else
			do i=1,p
				if ( mod(i,2) == 0) then
					call dgemm('t', 'n', n, l, m, 1.d0, a, m, uout, m, 0.d0, vout, n)
				else
					call dgemm('n', 'n', m, l, n, 1.d0, a, m, vout, n, 0.d0, uout, m)
				end if
			end do
		end if
		
		
		!Mise a jour de currentS !
		if (mod(p,2) == 0) then 
			if (currentS == 1) then
				currentS = 0
			else
				currentS = 1
			end if
		else
			currentS = currentS
		end if
		
		!Orthonormalisation!
		if (currentS ==0) then
			call orth_basis_ip(vout, n, n, l)
		else
			call orth_basis_ip(uout, m, m, l)

                end if
		
		if (currentS == 1) then
			! z=a'*u
			call dgemm('t', 'n', n, l, m, 1.d0, a, m, uout, m, 0.d0, z, n)
			!!compute  y=a*z
			call dgemm('n', 'n', m, l, n, 1.d0, a, m, z, n, 0.d0, y, m)
			!!compute  h=u'*y
		        call dgemm('t', 'n', l, l, m, 1.d0, uout, m, y, m, 0.d0, h, l)

		else
			! y=a*v
			call dgemm('n', 'n', m, l, n, 1.d0, a, m, vout, n, 0.d0, y, m)
			!!compute  z=a'*y
			call dgemm('t', 'n', n, l, m, 1.d0, a, m, y, m, 0.d0, z, n)
			!!compute  h=v'*z
			call dgemm('t', 'n', l, l, n, 1.d0, vout, n, z, n, 0.d0, h, l)
                       
		end if


		!décomposition spectrale!
		!Utilisation de z comme espace de stockage!
		call dsyev('v', 'u', l, h, l, w, z, lwork, ierr)

               
          ! Utilisation de la premiere colonne de z comme espace de stockage
           do i=1,(l)/2
          tmp        = w(i)
          w(i)       = w(l-i +1)
          w(l-i+1)   = tmp
          ! use the first column of z as a workspace
          z(1:l,1)   = h(:,i)
          h(:,i)     = h(:,l-i +1)
          h(:,l-i+1) = z(1:l,1)
       end do

		
		!V = VX 
		if (currentS == 1) then
			!Récupération de v à partir de u!
			call dgemm('t', 'n', n, l, m, 1.d0, a, m, uout, m, 0.d0, vout, n)
		else
			!Récupération de u à partir de v!
			call dgemm('n', 'n', m, l, n, 1.d0, a, m, vout, n, 0.d0, uout, m)
		end if
		!v=v*h
	        call dgemm('n', 'n', n, l, l, 1.d0, vout, n, h, l, 0d0, z, n)
                vout = z

		do i = (n_sv + 1),m
                 conv = 0
		  ! calcul de la norme!
		  !utilisation de Z comme espace de stockage!
		  if( i .gt. l) then
			exit
		 else
                        if (currentS == 1) then
                          normeA = (traceAM)
                        else
                          normeA = (traceAN)
                        end if
                         call dgemv('n', m, n, 1.d0, a, m, vout(:,i), 1, 0.d0, y, 1)
                        z(:,1) = vout(:,i)
                        call dgemv('t', m, n, 1.d0, a, m, y, 1, -w(i), z, 1)   
                        res = dnrm2(n, z, 1)
                         
			res = res/normeA
			if (res.gt.eps) then
				exit
			else
			        conv= conv + 1
			        res_sv(i)    = res
			        it_sv(i)     = it
                percentReached = 1 - sqrt(w(i)/w(1))
			end if
		end if
          n_sv = n_sv + conv
       end do

	write(*,'(" IT:",i5," -- Found ",i4," singular values (",f8.5,"% achieved)",a)',advance='no')it,n_sv,percentReached,char(13)
	write(*,'(" ")')
	it = it + 1
	if (percentReached .gt. percentage) then 
          exit
        end if
        if (n_sv .ge. l) then 
          exit
        end if
       
	end do
		
	n_sv = max(n_sv-1,1)
	
	!Racine carré des vs!
	 do j = 1, n_sv
        s(j) = sqrt(w(j))
    end do
    
	!Récupération de Vout et Uout!
	if (currentS == 1) then
        u(:,1:n_sv)= uout(:,1:n_sv)
		call dgemm('t', 'n', n_sv, l, m, 1.d0, a, m, y, m, 0.d0, v, n)
		!Normalisation!
		do j = 1, n_sv
			v(1:n,j) = v(1:n,j)/s(j)
		end do
	else
        v(:,1:n_sv) = vout(:,1:n_sv)
		call dgemm('n', 'n', m, l, n_sv, 1.d0, a, m, z, n, 0.d0, u, m)
        v = z
		!Normalisation!
		do j = 1, n_sv
			u(1:m,j) = u(1:m,j)/s(j)
		end do
	end if
        
	deallocate(uout,vout,z,h,w,y,matM, matN)
    return
    

  end subroutine subspace_iter_sv_both_alternate

  



  !!========================================================================
  !! orthogonalization using gram schmidt procedure
  !!========================================================================
  subroutine orth_basis(u_in,n,m,u_out)
    implicit none
    integer          ,intent(in)                            :: n,m
    double precision,dimension(n,m),intent(in)              :: u_in
    double precision,dimension(n,m),intent(out)             :: u_out
    !! local variables
    integer ::i,j
    do i=1,m
       u_out(:,i)=u_in(:,i)
    end do
    do j=1,m
       if( j.gt. 1)then
          do i=1,j-1
             u_out(:,j)= u_out(:,j) - dot_product(u_out(:,j),              &
                  u_out(:,i))*u_out(:,i)
          end do
       end if
       u_out(:,j)=u_out(:,j)/sqrt(dot_product(u_out(:,j),u_out(:,j)))
    end do
  end subroutine orth_basis




  !!========================================================================
  !! orthogonalization using gram schmidt procedure
  !!========================================================================
  subroutine orth_basis_ip(u, ldu, m, n)
    implicit none
    integer                    :: m, n, ldu
    double precision           :: u(ldu,n)

    !! local variables
    double precision, external :: dnrm2
    integer, parameter         :: ione=1
    integer                    :: k, j
    do k=1,n
       if( k.gt. 1)then
          do j=1,k-1
             u(1:m,k)= u(1:m,k) - u(1:m,j)*dot_product(u(1:m,k), u(1:m,j))
          end do
       end if
       u(1:m,k)=u(1:m,k)/dnrm2(m, u(1,k), ione)
    end do
  end subroutine orth_basis_ip

  






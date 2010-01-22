module utilities

  real*8 :: pi
  parameter (pi = 4.0d0*datan(1.0d0))

contains

  subroutine linear_solver(B, A, X, n)
    !
    !  Solves the linear system A X = B
    !
    implicit none
    
    integer, intent(in) :: n
    real*8, dimension(n), intent(inout) :: B
    real*8, dimension(n,n), intent(inout) :: A
    real*8, dimension(n), intent(out) :: X

    integer, dimension(n) :: pvt

    call legs(A, n, B, X, pvt)

  end subroutine linear_solver

  subroutine legs(A, n, B, X, indx)
    !
    ! Subroutine to solve the equation A(n,n)*X(n) = B(n) with the
    ! partial-pivoting Gaussian elimination scheme.
    !
    implicit none

    integer, intent (in) :: n
    real*8, dimension (n), intent (inout) :: B
    real*8, dimension (n,n), intent (inout) :: A
    integer, dimension (n), intent (out) :: indx
    real*8, dimension (n), intent (out) :: X

    integer :: i,j

    call elgs(A, n, indx)

    do i=1,n-1
       do j=i+1,n
          B(indx(j)) = B(indx(j))-A(indx(j),i)*B(indx(i))
       end do
    end do
    
    X(n) = B(indx(n))/A(indx(n),n)
    do i=n-1,1,-1
       X(i) = B(indx(i))
       do j=i+1,n
          X(i) = X(i)-A(indx(i),j)*X(j)
       end do
       X(i) =  X(i)/A(indx(i),i)
    end do
    
  end subroutine legs
  
  subroutine elgs(A, n, indx)
    !
    ! Subroutine to perform the partial-pivoting Gaussian elimination.
    ! A(n,n) is the original matrix in the input and transformed matrix
    ! plus the pivoting element ratios below the diagonal in the output.
    ! indx(n) records the pivoting order.
    !
    implicit none

    integer, intent (in) :: n
    real*8, dimension (n,n), intent (inout) :: A
    integer, dimension (n), intent (out) :: indx

    integer :: i,j,k,itmp
    real*8 :: C1,pi0,pi1,pj
    real*8, dimension (n) :: C

    !
    ! Initialize the index
    !
    do i=1,n
       indx(i) = i
    end do
    !
    ! Find the rescaling factors, one from each row
    !
    do i=1,n
       C1= 0.0
       do j=1,n
          C1 = amax1(C1,abs(A(i,j)))
       end do
       C(i) = C1
    end do
    !
    ! Search the pivoting (largest) element from each column
    !
    do j=1,n-1
       pi1 = 0.0
       do i=j,n
          pi0 = abs(A(indx(i),j))/C(indx(i))
          if(pi0.gt.pi1) then
             pi1 = pi0
             k   = i
          endif
       end do
       !
       ! Interchange the rows via INDX(N) to record pivoting order
       !
       itmp    = indx(j)
       indx(j) = indx(k)
       indx(k) = itmp
       do i=j+1,n
          pj  = A(indx(i),j)/A(indx(j),j)
          !
          ! Record pivoting ratios below the diagonal
          !
          A(indx(i),j) = pj
          !
          ! Modify other elements accordingly
          !
          do k=j+1,n
             A(indx(i),k) = A(indx(i),k)-pj*A(indx(j),k)
          end do
       end do
    end do
    
  end subroutine elgs

  subroutine interpolant(x,y,xi,yi)
    !
    !
    !
    implicit none

    real *8, intent(in) :: xi
    real *8, dimension(:), intent(in) :: x,y
    real *8, intent(out) :: yi

    integer :: nx,i

    nx = size(x)

    if(xi.lt.x(1)) then
       yi = y(1)
       write(*,*) ' WARNING!!! Interpolation on array out of scope 1'
       return
    end if

    if(xi.gt.x(nx)) then
       yi = y(nx)
       write(*,*) ' WARNING!!! Interpolation on array out of scope 2'
       return
    end if

    i = iminloc(dabs(x-xi))

    if(x(i).eq.xi) then
       yi = y(i)
       return
    end if

    if(x(i).gt.xi) i = i-1

    yi = y(i)+(xi-x(i))*(y(i+1)-y(i))/(x(i+1)-x(i))

  end subroutine interpolant

  function iminloc(arr)
    !
    !
    !
    implicit none

    real*8, dimension(:), intent(in) :: arr
    integer :: iminloc

    integer, dimension(1) :: imin

    imin = minloc(arr(:))

    iminloc = imin(1)

  end function iminloc

  function normvec(vec, n, p)
    !
    !  Computes the norm-p of a vector
    !
    implicit none

    integer, intent(in) :: n,p
    real*8, dimension(n) :: vec
    real*8 :: normvec

    integer :: i

    forall(i=1:n)
       vec(i) = vec(i)**p
    end forall

    normvec = sum(vec)
    normvec = normvec**(1./p)

  end function normvec

end module utilities

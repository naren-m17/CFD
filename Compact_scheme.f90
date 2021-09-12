! ------------------------------------------------------------------------------------------!
! 					                Name : Naren											!
!    						M. Tech in Computational Mechanics								!
!    						Department of Mechanical Engineering							!
!    								IIT- Guwahati											!
!																							!
!																							!
!    Compact Difference Scheme				    											!
!																							!
! ------------------------------------------------------------------------------------------!

!to solve the eq Ax = b
module solver
    implicit none
    
contains
!to eliminate the first elements from second row
!input is taken as 3 x n matrix A and 1 x n b vector
subroutine elimination(A, b)
implicit none
real(kind=8), DIMENSION(:,:), INTENT(INOUT) :: A !define all variables
real(kind=8), DIMENSION(:,:), INTENT(INOUT) :: b
INTEGER,dimension(2) :: sizeb
integer ::i
real (kind=8):: fact

sizeb = shape(b) !get the shape of b matrix

do i=2,sizeb(1)
    fact=A(i,1)/A(i-1,2)
    A(i,2)=A(i,2)-fact*A(i-1,3) !substract the lower diagonal element
    b(i,:)=b(i,:)-fact*b(i-1,:) 
end do

end subroutine elimination

!to eliminate the first elements from second row (pentadigonal matrix)
!input is taken as n x 5 matrix A and 1 x n b vector
subroutine elimination5(A, b)
    implicit none
    real(kind=8), DIMENSION(:,:), INTENT(INOUT) :: A !declare the variables
    real(kind=8), DIMENSION(:,:), INTENT(INOUT) :: b
    INTEGER,dimension(2) :: sizeb
    integer ::i
    real (kind=8):: fact
    
    sizeb = shape(b) !get the shape

    
    fact=A(2,2)/A(1,3) !eliminate the A(2,1)
    A(2,2:4)=A(2,2:4)-fact*A(1,3:5)
    
    do i=3,sizeb(1)
        fact=A(i,1)/A(i-2,3) !eliminate the A(i,i-2)
        A(i,1:3)=A(i,1:3)-fact*A(i-2,3:5)
        b(i,:)=b(i,:)-fact*b(i-2,:)

        fact=A(i,2)/A(i-1,3) !eliminate the A(i,i-1)
        A(i,2:4)=A(i,2:4)-fact*A(i-1,3:5)
        b(i,:)=b(i,:)-fact*b(i-1,:)
    end do
end subroutine elimination5    

!to back substitute and get the result
function backSub(A,b) result(x)
    implicit none
    real(kind=8), DIMENSION(:,:), INTENT(in) :: A !declare variable
    real(kind=8), DIMENSION(:,:), INTENT(in) :: b
    real(kind=8), DIMENSION(:,:), allocatable :: x
    INTEGER,dimension(2):: sizeb
    INTEGER ::i

    sizeb=shape(b) !get the shape of b
    allocate(x(sizeb(1),sizeb(2))) !allocate x matrix same as b

    x(sizeb(1),:)=b(sizeb(1),:)/A(sizeb(1),2) !find the solution of last 

    do i=sizeb(1)-1,1,-1
        x(i,:)=(b(i,:)-A(i,3)*x(i+1,:))/A(i,2) !find the all other solution
    end do

end function backSub


!to back substitute and get the result
function backSub5(A,b) result(x)
    implicit none
    real(kind=8), DIMENSION(:,:), INTENT(in) :: A !declare the variables
    real(kind=8), DIMENSION(:,:), INTENT(in) :: b
    real(kind=8), DIMENSION(:,:), allocatable :: x
    INTEGER,dimension(2):: sizeb
    INTEGER :: i

    sizeb=shape(b) !get the shape of b matrix
    
    allocate(x(sizeb(1),sizeb(2))) !allocate the x same as b

    x(sizeb(1),:)=b(sizeb(1),:)/A(sizeb(1),3) !find the solution of last element
    x(sizeb(1)-1,:)=(b(sizeb(1),:)-A(sizeb(1)-1,4)*x(sizeb(1),:))/A(sizeb(1),3) !find second to last x

    do i=sizeb(1)-2,1,-1
        x(i,:)=(b(i,:)-A(i,4)*x(i+1,:)-A(i,5)*x(i+2,:))/A(i,3) !find all other x values
    end do

end function backSub5

!to do TDMA and solve
function tdma(A,b) result(x)
    implicit none
    real(kind=8), DIMENSION(:,:), INTENT(inout) :: A !declare the varialbes
    real(kind=8), DIMENSION(:,:), INTENT(inout) :: b
    real(kind=8), DIMENSION(:,:), allocatable :: x
    
    call elimination(A,b) !elimination lower triangular matrix
    
    x=backSub(A,b) !do back substitution

end function tdma


!to do pDMA and solve
function pdma(A,b) result(x)
    implicit none
    real(kind=8), DIMENSION(:,:), INTENT(inout) :: A 
    real(kind=8), DIMENSION(:,:), INTENT(inout) :: b
    real(kind=8), DIMENSION(:,:), allocatable :: x
    
    call elimination5(A,b)
    
    x=backSub5(A,b)

end function pdma

function gaussSeidel(A,b,tol) result(x)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: A !declare all variables
    real(kind=8), dimension(:,:), intent(in) :: b
    real(kind=8), dimension(:,:),allocatable :: x
    real(kind=8), intent(in) :: tol
    real(kind=8),dimension(:), allocatable :: sum1
    INTEGER,dimension(2):: sizeb
    integer :: i,j

    sizeb=shape(b)!get shape

    allocate(x(sizeb(1),sizeb(2))) !allocate x
    allocate(sum1(sizeb(2))) !allocate sum
    x= 0!initial values

    do while(sum(abs(matmul(A,x)-b))/(sizeb(1)**2)>tol)
        do i=1,sizeb(1)
            sum1=0.0 !set to zero
            do j=1,sizeb(1)
                if(i.ne.j) then
                    sum1=sum1-A(i,j)*x(j,:) !add elements that are not zero
                end if
            end do
            x(i,:)=(sum1+b(i,:))/A(i,i) !find the solution using the formula
        end do
    end do

end function gaussSeidel
    
end module solver

module compactDiff
    use solver
    implicit none
    
contains

!function to calculate fourth order derivatives and write to file (both d/dx and d/dy)
subroutine compDiff4(phi,  alpha, deltaX, deltaY)
    implicit none
    real(kind=8),dimension(:,:),intent(inout) :: phi !declare the variable
    real(kind=8),intent(in) :: alpha, deltaX,deltaY
    real(kind=8),dimension(:,:),allocatable :: ddx,ddy, A,b, temp
    real(kind=8) :: rhs1,rhs2,xLen,yLen
    INTEGER,dimension(2) :: sizePhi
    integer :: i,j
    
    sizePhi = shape(phi) !get the shape
    allocate(ddx(sizePhi(1),sizePhi(2))) !allocate ddx and ddy
    allocate(ddy(sizePhi(1),sizePhi(2)))
    xLen=deltaX*(sizePhi(1)-1) !find length of domain (not passed through the function)
    yLen=deltaY*(sizePhi(2)-1)
    
    rhs1 =2.0/3.0*(alpha+2) !find rhs side variable
    rhs2=1.0/3.0*(4*alpha-1)
    
    !calculate d/dx
    allocate(b(sizePhi(1)-4,sizePhi(2)))
    b(1:,1:)=rhs1/(2*deltaX)*(phi(4:sizePhi(1)-1,:)-phi(2:sizePhi(1)-3,:))+&
                rhs2/(4*deltaX)*(phi(5:sizePhi(1),:)-phi(1:sizePhi(1)-4,:)) !find rhs side of equation 
    
    allocate(A(sizePhi(1)-4,3)) !Find the A matrix in n x 3
    A(:, 1)=alpha
    A(:,2)=1
    A(:,3)=alpha
    
    !forward difference to i=1,2
    ddx(1:2,:)=1.0/(12.0*deltaX)*(-25.0*phi(1:2,:)+48.0*phi(2:3,:)-36.0*phi(3:4,:)+16*phi(4:5,:)-3*phi(5:6,:))
    !backward difference to i=m,m-1
    ddx(sizePhi(1)-1:sizePhi(1),:)=1.0/(12.0*deltaX)*(25.0*phi(sizePhi(1)-1:sizePhi(1),:)-&
    48.0*phi(sizePhi(1)-2:sizePhi(1)-1,:)+36.0*phi(sizePhi(1)-3:sizePhi(1)-2,:)-&
    16.0*phi(sizePhi(1)-4:sizePhi(1)-5,:)+3.0*phi(sizePhi(1)-5:sizePhi(1)-6,:))
    
    !substract the boundary condition
    b(1,:)=b(1,:)-alpha*ddx(2,:)
    b(sizePhi(1)-4,:)=b(sizePhi(1)-4,:)-alpha*ddx(sizePhi(1)-1,:)
    
    ddx(3:sizePhi(1)-2,:)=tdma(A,b)!find the solution for all variable
    
    OPEN(1, FILE = "dPhidx.dat", STATUS='unknown') !write to file for contour
    
    do i=1,sizePhi(1)
        do j=1,sizePhi(2)
            write(1,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
        end do
        write(1,*) 
    end do
    close(1)

    OPEN(3, FILE = "dPhidx4x5.dat", STATUS='unknown') !write to file for x=5

    i=(sizePhi(1)+1)/2
    do j=1,sizePhi(2)
        write(3,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
    end do
    write(3,*) 

    close(3)

    OPEN(4, FILE = "dPhidx4y5.dat", STATUS='unknown') !write to file for y=5

    j=(sizePhi(1)+1)/2
    do i=1,sizePhi(2)
        write(4,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
    end do
    write(4,*) 

    close(4)

    deallocate(b) !deallocate all used variable
    deallocate(A)
    deallocate(ddx)
    !calculate d/dy
    allocate(b(sizePhi(1),sizePhi(2)-4))
    b(1:,1:)=rhs1/(2*deltaY)*(phi(:,4:sizePhi(2)-1)-phi(:,2:sizePhi(2)-3))+&
                rhs2/(2*deltaY)*(phi(:,5:sizePhi(2))-phi(:,1:sizePhi(2)-4))
    
    allocate(A(sizePhi(2)-4,3))
    A(:, 1)=alpha
    A(:,2)=1
    A(:,3)=alpha
    
    ddy(:,1:2)=1.0/(12.0*deltaY)*(-25.0*phi(:,1:2)+48.0*phi(:,2:3)-36.0*phi(:,3:4)+16*phi(:,4:5)-3*phi(:,5:6))
    ddy(:,sizePhi(2)-1:sizePhi(2))=1.0/(12.0*deltaY)*(25.0*phi(:,sizePhi(2)-1:sizePhi(2))-&
    48.0*phi(:,sizePhi(2)-2:sizePhi(2)-1)+36.0*phi(:,sizePhi(2)-3:sizePhi(2)-2)-&
    16.0*phi(:,sizePhi(2)-4:sizePhi(2)-5)+3.0*phi(:,sizePhi(2)-5:sizePhi(2)-6))
    b(:,1)=b(:,1)-alpha*ddy(:,2)
    b(:,sizePhi(2)-4)=b(:,sizePhi(2)-4)-alpha*ddy(:,sizePhi(2)-1)
    allocate(temp(sizePhi(2)-4,sizePhi(1)))
    temp = transpose(b)
    
    ddy(:,3:sizePhi(2)-2)=transpose(tdma(A,temp))
    OPEN(2, FILE = "dPhidy.dat", STATUS='unknown')
    
    do i=1,sizePhi(1)
        do j=1,sizePhi(2)
            write(2,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddy(i,j)
        end do
        write(2,*) 
    end do
    close(2)
    
end subroutine compDiff4
    
!function to calculate sixth order derivatives and write to file (both d/dx and d/dy)
subroutine compDiff6(phi,  alpha, deltaX, deltaY)
    implicit none
    real(kind=8),dimension(:,:),intent(inout) :: phi !declare all variable
    real(kind=8),intent(in) :: alpha, deltaX,deltaY
    real(kind=8),dimension(:,:),allocatable :: ddx,ddy, A,b, temp
    real(kind=8) :: rhs1,rhs2,rhs3,xLen,yLen
    INTEGER,dimension(2) :: sizePhi
    integer :: i,j
        
    sizePhi = shape(phi) !get the shape of phi matrix
    allocate(ddx(sizePhi(1),sizePhi(2))) !allocate ddx and ddy same size as phi
    allocate(ddy(sizePhi(1),sizePhi(2)))
    xLen=deltaX*(sizePhi(1)-1)
    yLen=deltaY*(sizePhi(2)-1)
        
    rhs1 =1.0/6.0*(alpha+9) !find rhs coefficients
    rhs2=1.0/15.0*(32*alpha-9)
    rhs3=1.0/10.0*(-3*alpha+1)
        
    !calculate d/dx
    allocate(b(sizePhi(1)-6,sizePhi(2))) !allocate rhs side of the equation 
    b(1:,1:)=rhs1/(2*deltaX)*(phi(5:sizePhi(1)-2,:)-phi(3:sizePhi(1)-4,:))+&
                rhs2/(4*deltaX)*(phi(6:sizePhi(1)-1,:)-phi(2:sizePhi(1)-5,:))+&
                rhs3/(6*deltaX)*(phi(7:sizePhi(1),:)-phi(1:sizePhi(1)-6,:))
        
    allocate(A(sizePhi(1)-6,3)) !fill the A matrix
    A(:, 1)=alpha
    A(:,2)=1
    A(:,3)=alpha
        
    !forward difference and backward difference for first three and last three varibles
    ddx(1:3,:)=1.0/(12.0*deltaX)*(-25.0*phi(1:3,:)+48.0*phi(2:4,:)-36.0*phi(3:5,:)+16*phi(4:6,:)-3*phi(5:7,:))
    ddx(sizePhi(1)-2:sizePhi(1),:)=1.0/(12.0*deltaX)*(25.0*phi(sizePhi(1)-2:sizePhi(1),:)-&
        48.0*phi(sizePhi(1)-3:sizePhi(1)-1,:)+36.0*phi(sizePhi(1)-4:sizePhi(1)-2,:)-&
        16.0*phi(sizePhi(1)-5:sizePhi(1)-3,:)+3.0*phi(sizePhi(1)-6:sizePhi(1)-4,:))
        
    !change first and last elements 
    b(1,:)=b(1,:)-alpha*ddx(3,:)
    b(sizePhi(1)-6,:)=b(sizePhi(1)-6,:)-alpha*ddx(sizePhi(1)-2,:)
        
    ddx(4:sizePhi(1)-3,:)=tdma(A,b) !solve all other values of x
        
    OPEN(1, FILE = "dPhidx6.dat", STATUS='unknown') !write to file for ploting contour
        
    do i=1,sizePhi(1)
        do j=1,sizePhi(2)
            write(1,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
        end do
        write(1,*) 
    end do
    close(1)

    OPEN(3, FILE = "dPhidx6x5.dat", STATUS='unknown') !write to file at x=5

    i=(sizePhi(1)+1)/2
    do j=1,sizePhi(2)
        write(3,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
    end do
    write(3,*) 

    close(3)

    OPEN(4, FILE = "dPhidx6y5.dat", STATUS='unknown') !write to file at y=5

    j=(sizePhi(1)+1)/2
    do i=1,sizePhi(2)
        write(4,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
    end do
    write(4,*) 

    close(4)
        
    deallocate(b)
    deallocate(A)
    deallocate(ddx)
    !calculate d/dy
    allocate(b(sizePhi(1),sizePhi(2)-6))
    b(1:,1:)=rhs1/(2*deltaY)*(phi(:,5:sizePhi(2)-2)-phi(:,3:sizePhi(2)-4))+&
                rhs2/(4*deltaY)*(phi(:,6:sizePhi(2)-1)-phi(:,2:sizePhi(2)-5))+&
                rhs3/(6*deltaY)*(phi(:,7:sizePhi(2))-phi(:,1:sizePhi(2)-6))
    
    allocate(A(sizePhi(2)-6,3))
    A(:, 1)=alpha
    A(:,2)=1
    A(:,3)=alpha
        
    ddy(:,1:3)=1.0/(12.0*deltaY)*(-25.0*phi(:,1:3)+48.0*phi(:,2:4)-36.0*phi(:,3:5)+16*phi(:,4:6)-3*phi(:,5:7))
    ddy(:,sizePhi(2)-2:sizePhi(2))=1.0/(12.0*deltaY)*(25.0*phi(:,sizePhi(2)-2:sizePhi(2))-&
        48.0*phi(:,sizePhi(2)-3:sizePhi(2)-1)+36.0*phi(:,sizePhi(2)-4:sizePhi(2)-2)-&
        16.0*phi(:,sizePhi(2)-5:sizePhi(2)-3)+3.0*phi(:,sizePhi(2)-6:sizePhi(2)-4))

    b(:,1)=b(:,1)-alpha*ddy(:,3)
    b(:,sizePhi(2)-6)=b(:,sizePhi(2)-6)-alpha*ddy(:,sizePhi(2)-2)
    allocate(temp(sizePhi(2)-6,sizePhi(1)))
    temp = transpose(b)
        
    ddy(:,4:sizePhi(2)-3)=transpose(tdma(A,temp))
    OPEN(2, FILE = "dPhidy6.dat", STATUS='unknown')
        
    do i=1,sizePhi(1)
        do j=1,sizePhi(2)
            write(2,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddy(i,j)
        end do
        write(2,*) 
    end do
    close(2)
        
end subroutine compDiff6
    

!function to calculate eighth order derivatives and write to file (both d/dx and d/dy) using Gausss Seidel method
subroutine compDiff8G(phi,  alpha, deltaX, deltaY,tol)
    implicit none
    real(kind=8),dimension(:,:),intent(inout) :: phi !declare all variables
    real(kind=8),intent(in) :: alpha, deltaX,deltaY,tol
    real(kind=8),dimension(:,:),allocatable :: ddx, A,b
    real(kind=8) :: rhs1,rhs2,rhs3,xLen,yLen,beta
    INTEGER,dimension(2) :: sizePhi
    integer :: i,j
            
    sizePhi = shape(phi)
    allocate(ddx(sizePhi(1),sizePhi(2))) !allocate ddx same size as phi
    xLen=deltaX*(sizePhi(1)-1)
    yLen=deltaY*(sizePhi(2)-1)
            
    beta =1.0/20*(-3+alpha*8) !beta and rhs coefficients
    rhs1 =1.0/6.0*(12-7*alpha)
    rhs2=1.0/150.0*(568*alpha-183)
    rhs3=1.0/50.0*(9*alpha-4)
            
    !calculate d/dx
    allocate(b(sizePhi(1)-6,sizePhi(2))) !rhs side of the equation 
    b(1:,1:)=rhs1/(2*deltaX)*(phi(5:sizePhi(1)-2,:)-phi(3:sizePhi(1)-4,:))+&
                rhs2/(4*deltaX)*(phi(6:sizePhi(1)-1,:)-phi(2:sizePhi(1)-5,:))+&
                rhs3/(6*deltaX)*(phi(7:sizePhi(1),:)-phi(1:sizePhi(1)-6,:))
            
    allocate(A(sizePhi(1)-6,sizePhi(1)-6)) !set A matrix
    do i=1,sizePhi(1)-6
        A(i,i)=1
        if(i-1>0) then
            A(i,i-1)=alpha
        end if
        if(i-2>0) then
            A(i,i-2)=beta
        end if
        if(i<sizePhi(1)-6)then
            A(i,i+1)=alpha
        end if
        if(i<sizePhi(1)-5) then
            A(i,i+2)=beta
        end if
    end do

    !find first three and last three elements using forward differenc and backward difference        
    ddx(1:3,:)=1.0/(12.0*deltaX)*(-25.0*phi(1:3,:)+48.0*phi(2:4,:)-36.0*phi(3:5,:)+16*phi(4:6,:)-3*phi(5:7,:))
    ddx(sizePhi(1)-2:sizePhi(1),:)=1.0/(12.0*deltaX)*(25.0*phi(sizePhi(1)-2:sizePhi(1),:)-&
            48.0*phi(sizePhi(1)-3:sizePhi(1)-1,:)+36.0*phi(sizePhi(1)-4:sizePhi(1)-2,:)-&
            16.0*phi(sizePhi(1)-5:sizePhi(1)-3,:)+3.0*phi(sizePhi(1)-6:sizePhi(1)-4,:))
            
    !take care of elements in first two and last two that are already found
    b(1,:)=b(1,:)-alpha*ddx(3,:)-beta*ddx(2,:)
    b(2,:)=b(2,:)-beta*ddx(3,:)
    b(sizePhi(1)-6,:)=b(sizePhi(1)-6,:)-alpha*ddx(sizePhi(1)-2,:)-beta*ddx(sizePhi(1)-1,:)
    b(sizePhi(1)-7,:)=b(sizePhi(1)-7,:)-beta*ddx(sizePhi(1)-2,:)

    ddx(4:sizePhi(1)-3,:)=gaussSeidel(A,b,tol) !find the solution
        
    OPEN(1, FILE = "dPhidx8.dat", STATUS='unknown') !write to file for countour
            
    do i=1,sizePhi(1)
        do j=1,sizePhi(2)
            write(1,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
        end do
        write(1,*) 
    end do
    close(1)
    
    OPEN(3, FILE = "dPhidx8x5.dat", STATUS='unknown') !write to file at x=5

    i=(sizePhi(1)+1)/2
    do j=1,sizePhi(2)
        write(3,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
    end do
    write(3,*) 

    close(3)

    OPEN(4, FILE = "dPhidx8y5.dat", STATUS='unknown') !write to file at y=5

    j=(sizePhi(1)+1)/2
    do i=1,sizePhi(2)
        write(4,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
    end do
    write(4,*) 

    close(4)
            
    deallocate(b)
    deallocate(A)
    deallocate(ddx)
end subroutine compDiff8G

!function to calculate eighth order derivatives and write to file (both d/dx and d/dy) using Gausss Seidel method
subroutine compDiff8p(phi,  alpha, deltaX, deltaY)
    implicit none
    real(kind=8),dimension(:,:),intent(inout) :: phi 
    real(kind=8),intent(in) :: alpha, deltaX,deltaY
    real(kind=8),dimension(:,:),allocatable :: ddx,A,b
    real(kind=8) :: rhs1,rhs2,rhs3,xLen,yLen,beta
    INTEGER,dimension(2) :: sizePhi
    integer :: i,j
            
    sizePhi = shape(phi)
    allocate(ddx(sizePhi(1),sizePhi(2)))
    xLen=deltaX*(sizePhi(1)-1)
    yLen=deltaY*(sizePhi(2)-1)
            
    beta =1.0/20*(-3+alpha*8)
    rhs1 =1.0/6.0*(12-7*alpha)
    rhs2=1.0/150.0*(568*alpha-183)
    rhs3=1.0/50.0*(9*alpha-4)
            
    !calculate d/dx
    allocate(b(sizePhi(1)-6,sizePhi(2)))
    b=rhs1/(2*deltaX)*(phi(5:sizePhi(1)-2,:)-phi(3:sizePhi(1)-4,:))+&
                rhs2/(4*deltaX)*(phi(6:sizePhi(1)-1,:)-phi(2:sizePhi(1)-5,:))+&
                rhs3/(6*deltaX)*(phi(7:sizePhi(1),:)-phi(1:sizePhi(1)-6,:))
            
    allocate(A(sizePhi(1)-6,5))
    A(:,1)=beta
    A(:,2)=alpha
    A(:,3)=1
    A(:,4)=alpha
    A(:,5)=beta
        
    ddx(1:3,:)=1.0/(12.0*deltaX)*(-25.0*phi(1:3,:)+48.0*phi(2:4,:)-36.0*phi(3:5,:)+16*phi(4:6,:)-3*phi(5:7,:))
    ddx(sizePhi(1)-2:sizePhi(1),:)=1.0/(12.0*deltaX)*(25.0*phi(sizePhi(1)-2:sizePhi(1),:)-&
            48.0*phi(sizePhi(1)-3:sizePhi(1)-1,:)+36.0*phi(sizePhi(1)-4:sizePhi(1)-2,:)-&
            16.0*phi(sizePhi(1)-5:sizePhi(1)-3,:)+3.0*phi(sizePhi(1)-6:sizePhi(1)-4,:))
            
    b(1,:)=b(1,:)-alpha*ddx(3,:)-beta*ddx(2,:)
    b(2,:)=b(2,:)-beta*ddx(3,:)
    b(sizePhi(1)-6,:)=b(sizePhi(1)-6,:)-alpha*ddx(sizePhi(1)-2,:)-beta*ddx(sizePhi(1)-1,:)
    b(sizePhi(1)-7,:)=b(sizePhi(1)-7,:)-beta*ddx(sizePhi(1)-2,:)

    ddx(4:sizePhi(1)-3,:)=pdma(A,b)

    deallocate(b)
    deallocate(A)

    OPEN(1, FILE = "dPhidx8.dat", STATUS='unknown')
            
    do i=1,sizePhi(1)
        do j=1,sizePhi(2)
            write(1,*) xLen/(sizePhi(1)-1)*(i-1), yLen/(sizePhi(2)-1)*(j-1), ddx(i,j)
        end do
        write(1,*) 
    end do
    close(1)

    deallocate(ddx)
end subroutine compDiff8p

end module compactDiff

program main
    use compactDiff

    implicit none

    real (kind=8), DIMENSION(:,:), allocatable ::  phi !declare all variables
    real(kind=8) :: xLen = 10, yLen = 10, deltaX, deltaY
    INTEGER :: i, j, xGrid = 100, yGrid = 100, order

    write(*,*) "Select required order of accuracy (4,6,8): " !get the user input 
    read(*,*) order

    deltaX = xLen/xGrid
    deltaY = yLen/yGrid
    allocate(phi(xGrid+1,yGrid+1))

    do i = 1,xGrid+1
        do j = 1,yGrid+1
            phi(i,j)=sin(xLen/xGrid*(i-1))*cos(yLen/yGrid*(j-1)) !intial grid 
        end do
    end do

    OPEN(1, FILE = "InitialGrid.dat", STATUS='unknown') !contour the inital values

    do i=1,xGrid+1
        do j=1,yGrid+1
            write(1,*) xLen/xGrid*(i-1), yLen/yGrid*(j-1), phi(i,j)
        end do
        write(1,*) 
    end do

    close(1)

    !select and run the selected order
    if(order == 4) then
        call compDiff4(phi,dble(.3), deltaX, deltaY)
        write(*,*) "The derivatives are written to the file"
    else if(order==6) then
        call compDiff6(phi,dble(.5), deltaX, deltaY)
        write(*,*) "The derivatives are written to the file"
    else if(order==8) then
        call compDiff8G(phi,dble(.5), deltaX, deltaY,dble(1e-8))
        write(*,*) "The derivatives are written to the file"
    else 
        write(*,*) "Enter only (4/6/8) orders"
        write(*,*)
    end if

end program main
! ------------------------------------------------------------------------------------------!
! 					Name : Naren				Roll no.: 204103110							!
!    						M. Tech in Computational Mechanics								!
!    						Department of Mechanical Engineering							!
!    								IIT- Guwahati											!
!																							!
!    Programming assignment#2																!
!    Advanced Computational Fluid dynamics													!
!																							!
!    Solving Conduction equation (Laplace eq) using Incomplete LU   						!
!																							!
! ------------------------------------------------------------------------------------------!

subroutine exact()
    implicit none

    INTEGER::xGrid=60,yGrid=60,i,j,k
    double precision::xLen=1,yLen=1, pi, xDelta, yDelta, x, y, z

    xDelta=xLen/(xGrid-2)
    yDelta=yLen/(yGrid-2)
    pi=4.D0*DATAN(1.D0)

    OPEN(2,FILE="exact.dat", STATUS='unknown')
    OPEN(3,FILE="TexactX5.dat", STATUS='unknown')
    OPEN(4,FILE="TexactY5.dat", STATUS='unknown')

    do i=1,xGrid
        do j=1,yGrid
            x=xDelta*(i-1)
            y=yDelta*(j-1)
            z=0
            do k=1,110
                z=z+2.d0/pi*(((-1.d0)**(k+1)+1.d0)/k)*sin(k*pi*x)*sinh(k*pi*y)/sinh(k*pi)
                !write(*,*) z
            end do
            write(2,*)x, y, z
            if(i==xGrid/2 .and. modulo(j,4)==0 .and. j.ne.60)then
                write(3,*) y,z
            end if
            if(j==yGrid/2 .and. modulo(i,4)==0 .and. i.ne.60)then
                write(4,*) x,z
            end if
        end do
        write(2,*)
    end do

    CLOSE(2)
end subroutine exact

program ilu
    implicit none
    
    INTEGER :: xGrid=60,yGrid=60,i,j
    double precision :: xLen=1,yLen=1,xDelta,yDelta,beta,TL=0,TR=0,TB=0,TT=1, tolerance=1e-6,x,y
    double precision,DIMENSION(:,:), ALLOCATABLE::A,L,U,M
    double precision,DIMENSION(:),ALLOCATABLE::T,b, TForw, TNew

    call exact()

    xDelta=xLen/(xGrid-1)
    yDelta=yLen/(yGrid-1)
    beta=xDelta/yDelta

    allocate(A((xGrid-2)*(yGrid-2),5))
    allocate(L((xGrid-2)*(yGrid-2),3))
    allocate(U((xGrid-2)*(yGrid-2),2))
    allocate(M((xGrid-2)*(yGrid-2),7))
    allocate(T((xGrid-2)*(yGrid-2)))
    allocate(b((xGrid-2)*(yGrid-2)))
    allocate(TForw((xGrid-2)*(yGrid-2)))
    allocate(TNew((xGrid-2)*(yGrid-2)))

    A(:,1)=1
    A(:,2)=beta**2
    A(:,3)=(-2-2*beta**2)
    A(:,4)=beta**2
    A(:,5)=1
    L=0
    U=0
    T=0
    b=0
    TNew=1

    do i=1,xGrid-2
        do j=1,yGrid-2
            if(i==1) then
                A((i-1)*(xGrid-2)+j,1)=0
            else if(i==xGrid-2) then
                A((i-1)*(xGrid-2)+j,5)=0
            end if
            if(j==1) then
                A((i-1)*(xGrid-2)+j,2)=0
            else if(j==yGrid-2) then
                A((i-1)*(xGrid-2)+j,4)=0
            end if
        end do
    end do
    
    L(:,1)=A(:,1)
    L(:,2)=A(:,2)

    do i=1,xGrid-2
        do j=1,yGrid-2
            if(i==1 .and. j==1)then
                L(1,3)=A(1,3)
            else if(i==1) then
                L(j,3)=A(j,3)-L(j,2)*U(j-1,1)
            else
                L((i-1)*(xGrid-2)+j,3)=A((i-1)*(xGrid-2)+j,3)-L((i-1)*(xGrid-2)+j,1)*U((i-2)*(xGrid-2)+j,2)&
                                    -L((i-1)*(xGrid-2)+j,2)*U((i-1)*(xGrid-2)+j-1,1)
            end if
            U((i-1)*(xGrid-2)+j,1)=A((i-1)*(xGrid-2)+j,4)/L((i-1)*(xGrid-2)+j,3)
            U((i-1)*(xGrid-2)+j,2)=A((i-1)*(xGrid-2)+j,5)/L((i-1)*(xGrid-2)+j,3)
            if(j==1) then
                b((i-1)*(xGrid-2)+j)=-beta**2*TB
            else if(j==(yGrid-2)) then
                b((i-1)*(xGrid-2)+j)=-beta**2*TT
            end if
            if(i==xGrid-2)then
                b((i-1)*(xGrid-2)+j)=-TR
            end if
        end do
    end do

    b(1)=-TL-beta**2*TB
    b(2:(xGrid-2))=-TL
    b(xGrid-2)=-TL-beta**2*TT
    b((xGrid-2)*(yGrid-2))=-TR-beta**2*TT
    b((xGrid-3)*(yGrid-2)-1)=-TR-beta**2*TB

    L(1,2)=0
    L(1:xGrid-2,1)=0
    U((xGrid-2)*(yGrid-2),1)=0
    U((xGrid-3)*(yGrid-2)+1:(xGrid-2)*(yGrid-2),2)=0

    do i=1,(xGrid-2)*(yGrid-2)
        if(i==1)then
            M(i,3)=0
            M(i,4)=L(i,3)
            M(i,6)=0
        else if(i<=(xGrid-2)) then
            M(i,2)=0
            M(i,4)=L(i,2)*U(i-1,1)+L(i,3)
            M(i,6)=L(i,2)*U(i-1,2)
        else
            M(i,2)=L(i,1)*U(i-(xGrid-2),1)
            M(i,4)=L(i,1)*U(i-(xGrid-2),2)+L(i,2)*U(i-1,1)+L(i,3)
            M(i,6)=L(i,2)*U(i-1,2)
        end if

        M(i,1)=L(i,1)
        M(i,3)=L(i,2)
        M(i,5)=U(i,1)*L(i,3)
        M(i,7)=U(i,2)*L(i,3)
    end do

    open(1, FILE="Matrix.m",STATUS='unknown',recl=5000)

    write(1,*) "A=["
    do i=1,(xGrid-2)*(yGrid-2)
        write(1,*) A(i,:)
    end do
    write(1,*) "];"
    write(1,*) "L=["
    do i=1,(xGrid-2)*(yGrid-2)
        write(1,*) L(i,:)
    end do
    write(1,*) "];"
    write(1,*) "U=["
    do i=1,(xGrid-2)*(yGrid-2)
        write(1,*) U(i,:)
    end do
    write(1,*) "];"
    write(1,*) "M=["
    do i=1,(xGrid-2)*(yGrid-2)
        write(1,*) M(i,:)
    end do
    write(1,*) "];"
    write(1,*) "b=["
    do i=1,(xGrid-2)*(yGrid-2)
        write(1,*) b(i)
    end do
    write(1,*) "];"
    

    do while(sum(abs(T-TNew))/((xGrid-2)*(yGrid-2)) > tolerance)

    T=TNew

    TForw(1)=(b(1)+T(xGrid-2)*M(1,6))/L(1,3)
    do i=2,(xGrid-2)*(yGrid-2)
        if(i< xGrid-2) then
            TForw(i)=(b(i)+T(xGrid-3+i)*M(i,6)-L(i,2)*TForw(i-1))/L(i,3)
            !write(*,*) b(i)+T(xGrid-3+i)*M(i,6)
        else
            TForw(i)=(b(i)+T(xGrid-3+i)*M(i,6)+T(-xGrid+3+i)*M(i,2)-L(i,2)*TForw(i-1)-L(i,1)*TForw(i-xGrid+2))/L(i,3)
            !write(*,*)b(i)+T(xGrid-3+i)*M(i,6)+T(-xGrid+3+i)*M(i,2)
        end if
    end do

    !write(*,*) TForw(i)

    TNew((xGrid-2)*(yGrid-2))=TForw((xGrid-2)*(yGrid-2))
    do i=(xGrid-2)*(yGrid-2)-1,1,-1
        if(i > (xGrid-3)*(yGrid-2)) then
            TNew(i)=TForw(i)-U(i,1)*TNew(i+1)
        else
            TNew(i)=(TForw(i)-U(i,1)*TNew(i+1)-U(i,2)*TNew(i+xGrid-2))
        end if
    end do
    !write(*,*) TNew
    !exit
    !write(*,*) sum(abs(T-TNew))/((xGrid-2)*(yGrid-2))
end do
write(1,*) "T=["
do i=1,(xGrid-2)
    write(1,*) TNew((i-1)*(yGrid-2)+1:i*(yGrid-2))
end do
write(1,*) "];"
close(1)

OPEN(2,FILE="T.dat", STATUS='unknown')
OPEN(3,FILE="TX5.dat", STATUS='unknown')
OPEN(4,FILE="TY5.dat", STATUS='unknown')

do i=1,xGrid-2
    do j=1,yGrid-2
        x=xDelta*(i)
        y=yDelta*(j)
        write(2,*) x,y,TNew((i-1)*(xGrid-2)+j)
        if(i==xGrid/2)then
            write(3,*) y,TNew((i-1)*(xGrid-2)+j)
        end if
        if(j==yGrid/2)then
            write(4,*) x,TNew((i-1)*(xGrid-2)+j)
        end if
    end do
end do

CLOSE(2)

end program ilu
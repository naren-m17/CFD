! ------------------------------------------------------------------------------------------!
! 					Name : Naren				Roll no.: 204103110							!
!    						M. Tech in Computational Mechanics								!
!    						Department of Mechanical Engineering							!
!    								IIT- Guwahati											!
!																							!
!    Programming assignment#2																!
!    Advanced Computational Fluid dynamics													!
!																							!
!    Solving Conduction equation (Laplace eq) using Strongly Implicit Method   						!
!																							!
! ------------------------------------------------------------------------------------------!

subroutine exact()
    implicit none

    INTEGER::xGrid=100,yGrid=100,i,j
    double precision::xLen=2,yLen=2, pi, xDelta, yDelta, x, y, z, Tm=2.d0

    xDelta=xLen/(xGrid)
    yDelta=yLen/(yGrid)
    pi=4.D0*DATAN(1.D0)

    OPEN(2,FILE="exact.dat", STATUS='unknown')
    OPEN(3,FILE="TexactX5.dat", STATUS='unknown')
    OPEN(4,FILE="TexactY5.dat", STATUS='unknown')

    do i=1,xGrid
        do j=1,yGrid
            x=xDelta*(i-1)
            y=yDelta*(j-1)
            z=Tm*sin(pi*x/xLen)*sinh(pi*y/xLen)/sinh(pi*yLen/xLen)
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

program sip
    implicit none
    
    INTEGER :: xGrid=100,yGrid=100,i,j
    double precision :: xLen=2,yLen=2,TL=0,TR=0,TB=0,TT=2, tolerance=1e-6, alpha=0.5
    double precision :: x,xDelta,yDelta,beta, pi,y
    double precision,DIMENSION(:,:), ALLOCATABLE::A,L,U,M
    double precision,DIMENSION(:),ALLOCATABLE::T,b, TForw, TNew

    xDelta=xLen/(xGrid-1)
    yDelta=yLen/(yGrid-1)
    beta=xDelta/yDelta
    pi=4.D0*DATAN(1.D0)

    allocate(A((xGrid-2)*(yGrid-2),5))
    allocate(L((xGrid-2)*(yGrid-2),3))
    allocate(U((xGrid-2)*(yGrid-2),2))
    allocate(M((xGrid-2)*(yGrid-2),7))
    allocate(T((xGrid-2)*(yGrid-2)))
    allocate(b((xGrid-2)*(yGrid-2)))
    allocate(TForw((xGrid-2)*(yGrid-2)))
    allocate(TNew((xGrid-2)*(yGrid-2)))

    call exact()

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

    do i=1,xGrid-2
        do j=1,yGrid-2
            if(i==1 .and. j==1)then
                L((i-1)*(xGrid-2)+j,1)=A((i-1)*(xGrid-2)+j,1)
                L((i-1)*(xGrid-2)+j,2)=A((i-1)*(xGrid-2)+j,2)
                L((i-1)*(xGrid-2)+j,3)=A((i-1)*(xGrid-2)+j,3)
                U((i-1)*(xGrid-2)+j,1)=A((i-1)*(xGrid-2)+j,4)/L((i-1)*(xGrid-2)+j,3)
                U((i-1)*(xGrid-2)+j,2)=A((i-1)*(xGrid-2)+j,5)/L((i-1)*(xGrid-2)+j,3)
            else if(i==1) then
                L((i-1)*(xGrid-2)+j,1)=A((i-1)*(xGrid-2)+j,1)
                L((i-1)*(xGrid-2)+j,2)=A((i-1)*(xGrid-2)+j,2)/(1+alpha*U((i-1)*(xGrid-2)+j-1,2))
                L((i-1)*(xGrid-2)+j,3)=A((i-1)*(xGrid-2)+j,3)+alpha*(L((i-1)*(xGrid-2)+j,2)*U((i-1)*(xGrid-2)+j-1,2))-&
                                        L((i-1)*(xGrid-2)+j,2)*U((i-1)*(xGrid-2)+j-1,1)
                U((i-1)*(xGrid-2)+j,1)=A((i-1)*(xGrid-2)+j,4)/L((i-1)*(xGrid-2)+j,3)
                U((i-1)*(xGrid-2)+j,2)=(A((i-1)*(xGrid-2)+j,5)-alpha*L((i-1)*(xGrid-2)+j,2)*U((i-1)*(xGrid-2)+j-1,2))&
                                        /L((i-1)*(xGrid-2)+j,3)
            else
                L((i-1)*(xGrid-2)+j,1)=A((i-1)*(xGrid-2)+j,1)/(1+alpha*U((i-2)*(xGrid-2)+j,1))
                L((i-1)*(xGrid-2)+j,2)=A((i-1)*(xGrid-2)+j,2)/(1+alpha*U((i-1)*(xGrid-2)+j-1,2))
                L((i-1)*(xGrid-2)+j,3)=A((i-1)*(xGrid-2)+j,3)+alpha*(L((i-1)*(xGrid-2)+j,1)*U((i-2)*(xGrid-2)+j,1)+&
                                    L((i-1)*(xGrid-2)+j,2)*U((i-1)*(xGrid-2)+j-1,2))-&
                            L((i-1)*(xGrid-2)+j,1)*U((i-2)*(xGrid-2)+j,2)-L((i-1)*(xGrid-2)+j,2)*U((i-1)*(xGrid-2)+j-1,1)
                U((i-1)*(xGrid-2)+j,1)=(A((i-1)*(xGrid-2)+j,4)-alpha*L((i-1)*(xGrid-2)+j,1)*U((i-2)*(xGrid-2)+j,1))&
                                        /L((i-1)*(xGrid-2)+j,3)
                U((i-1)*(xGrid-2)+j,2)=(A((i-1)*(xGrid-2)+j,5)-alpha*L((i-1)*(xGrid-2)+j,2)*U((i-1)*(xGrid-2)+j-1,2))&
                                        /L((i-1)*(xGrid-2)+j,3)
            end if
            if(j==1) then
                b((i-1)*(xGrid-2)+j)=-beta**2*TB
            else if(j==(yGrid-2)) then
                x=xDelta*(i-1)
                b((i-1)*(xGrid-2)+j)=-beta**2*TT*sin(pi*x/xLen)
            end if
            if(i==xGrid-2)then
                b((i-1)*(xGrid-2)+j)=-TR
            end if
        end do
    end do

    b(1)=-TL-beta**2*TB
    b(2:(xGrid-2))=-TL
    b(xGrid-2)=-TL-beta**2*TT*sin(pi*0/xLen)
    b((xGrid-2)*(yGrid-2))=-TR-beta**2*TT*sin(pi*2.d0/xLen)
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

    do i=1,xGrid-2
        do j=1,yGrid-2
            if(j==1) then
                M((i-1)*(xGrid-2)+j,2)=0
            else if(j==xGrid-2) then
                M((i-1)*(xGrid-2)+j,6)=0
            end if
        end do
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
            TForw(i)=(b(i)+(T(xGrid-3+i)-alpha*(T(xGrid-2+i)+T(i+1)-T(i)))*M(i,6)-L(i,2)*TForw(i-1))/L(i,3)
            !write(*,*) b(i)+T(xGrid-3+i)*M(i,6)
        else
            TForw(i)=(b(i)+(T(xGrid-3+i)-alpha*(T(xGrid-2+i)+T(i+1)-T(i)))*M(i,6)+(T(-xGrid+3+i)-&
                    alpha*(T(-xGrid+2+i)+T(i-1)-T(i)))*M(i,2)-L(i,2)*TForw(i-1)-L(i,1)*TForw(i-xGrid+2))/L(i,3)
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
        if(i==(xGrid-2)/2)then
            write(3,*) y,TNew((i-1)*(xGrid-2)+j)
        end if
        if(j==(yGrid-2)/2)then
            write(4,*) x,TNew((i-1)*(xGrid-2)+j)
        end if
    end do
    write(2,*)
end do


end program sip
! ------------------------------------------------------------------------------------------!
! 					Name : Naren				Roll no.: 204103110							!
!    						M. Tech in Computational Mechanics								!
!    						Department of Mechanical Engineering							!
!    								IIT- Guwahati											!
!																							!
!    Programming assignment#2																!
!    Advanced Computational Fluid dynamics													!
!																							!
!    Solving Conduction equation (Laplace eq) using conjugate gradient method  				!
!																							!
! ------------------------------------------------------------------------------------------!

subroutine exact()
    implicit none

    INTEGER::xGrid=100,yGrid=100,i,j
    double precision::xLen=1,yLen=1, pi, xDelta, yDelta, x, y

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
            write(2,*)x, y, sin(2*pi*x)*sin(2*pi*y)
            if(i==xGrid/2 .and. modulo(j,4)==0 .and. j.ne.60)then
                write(3,*) y,sin(2*pi*x)*sin(2*pi*y)
            end if
            if(j==yGrid/2 .and. modulo(i,4)==0 .and. i.ne.60)then
                write(4,*) x,sin(2*pi*x)*sin(2*pi*y)
            end if
        end do
        write(2,*)
    end do

    CLOSE(2)
end subroutine exact

program conjugate
    implicit none
    
    INTEGER :: xGrid=80,yGrid=80,i,j
    double precision :: xLen=1,yLen=1,xDelta,yDelta,beta,beta1, alpha,tolerance=1e-6, pi,x,y
    double precision,DIMENSION(:,:), ALLOCATABLE::A
    double precision,DIMENSION(:),ALLOCATABLE::T,b, r, rp, p

    xDelta=xLen/(xGrid-1)
    yDelta=yLen/(yGrid-1)
    beta=xDelta/yDelta
    pi=4.D0*DATAN(1.D0)

    allocate(A((xGrid-2)*(yGrid-2),(xGrid-2)*(yGrid-2)))
    allocate(b((xGrid-2)*(yGrid-2)))
    allocate(r((xGrid-2)*(yGrid-2)))
    allocate(rp((xGrid-2)*(yGrid-2)))
    allocate(p((xGrid-2)*(yGrid-2)))
    allocate(T((xGrid-2)*(yGrid-2)))
    call exact()

    A(2,1)=beta**2
    A((xGrid-2)*(yGrid-2),(xGrid-2)*(yGrid-2)-1)=beta**2

    OPEN(1,file='Matrix.m',STATUS='unknown', recl=175000)

    do i=1,xGrid-2
        do j=1,yGrid-2
            if(i==1) then
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j)=(-2-2*beta**2)
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j+1)=beta**2
                A((i-1)*(xGrid-2)+j,i*(xGrid-2)+j)=1.d0
                if(j>1) then
                    A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j-1)=beta**2
                end if
                
            else if(i==xGrid-2) then
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j)=(-2-2*beta**2)            
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j-1)=beta**2
                A((i-1)*(xGrid-2)+j,(i-2)*(xGrid-2)+j)=1.d0
                if(j<xGrid-2) then
                    A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j+1)=beta**2
                end if
                
            else
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j)=(-2-2*beta**2)
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j-1)=beta**2
                A((i-1)*(xGrid-2)+j,(i-2)*(xGrid-2)+j)=1.d0
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j+1)=beta**2
                A((i-1)*(xGrid-2)+j,i*(xGrid-2)+j)=1.d0

            end if

            if(j==1 .and. i>1) then
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j-1)=0
            else if(j==yGrid-2 .and. i<xGrid-2) then
                A((i-1)*(xGrid-2)+j,(i-1)*(xGrid-2)+j+1)=0
            end if
            x=(xDelta-2)*(i-1)
            y=(yDelta-2)*(j-1)
            b((i-1)*(xGrid-2)+j)=-8*pi**2*sin(2*pi*x)*sin(2*pi*y)*xDelta**2
        end do
    end do
    
    write(1,*) "A=["
    do i=1,(xGrid-2)*(yGrid-2)
        !write(1,*) A(i,:)
    end do
    write(1,*) "];"

    write(1,*) "b=["
    do i=1,(xGrid-2)*(yGrid-2)
        write(1,*) b(i)
    end do
    write(1,*) "];"

    x=0
    r=b
    p=0
    i=0
    do while(maxval(abs(r))>tolerance)
       if(i==0) then
        beta1=0
       else
        beta1=dot_product(r,r)/dot_product(rp,rp)
       end if
       p=r+beta1*p
       alpha=dot_product(r,r)/dot_product(p,reshape(matmul(A,reshape(p,(/(xGrid-2)*(yGrid-2),1/))),(/(xGrid-2)*(yGrid-2)/)))
       T=T+alpha*p
       rp=r
       r=b-reshape(matmul(A,reshape(T,(/(xGrid-2)*(yGrid-2),1/))),(/(xGrid-2)*(yGrid-2)/))
       i=i+1
       write(*,*) maxval(abs(r)), beta1, alpha
    end do

    write(1,*) "T=["
    do i=1,(xGrid-2)
        write(1,*) T((i-1)*(yGrid-2)+1:i*(yGrid-2))
    end do
    write(1,*) "];"

    CLOSE(1)

    OPEN(2,FILE="T.dat", STATUS='unknown')
    OPEN(3,FILE="TX5.dat", STATUS='unknown')
    OPEN(4,FILE="TY5.dat", STATUS='unknown')

    do i=1,xGrid-2
        do j=1,yGrid-2
            x=xDelta*(i)
            y=yDelta*(j)
            write(2,*) x,y,T((i-1)*(xGrid-2)+j)
            if(i==(xGrid-2)/2)then
                write(3,*) y,T((i-1)*(xGrid-2)+j)
            end if
            if(j==(yGrid-2)/2)then
                write(4,*) x,T((i-1)*(xGrid-2)+j)
            end if
        end do
        write(2,*)
    end do

end program conjugate
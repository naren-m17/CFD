! ------------------------------------------------------------------------------------------!
! 					Name : Naren				Roll no.: 204103110							!
!    						M. Tech in Computational Mechanics								!
!    						Department of Mechanical Engineering							!
!    								IIT- Guwahati											!
!																							!
!    Programming assignment#1																!
!    Advanced Computational Fluid dynamics													!
!																							!
!    Flow in 2D lid-driven cavity using MAC													!
!																							!
! ------------------------------------------------------------------------------------------!

program MAC
    implicit none
    
    INTEGER   :: xGrid = 300, yGrid = 300, i, j
    double precision, DIMENSION(:,:), ALLOCATABLE :: u, v, p, uNew, vNew, pNew, u2, uv, v2, rhsu, rhsv, psi, um, vm
    double precision :: uVelTop = 1, tolerance = 1e-7, xlen = 1.0, ylen = 1.0, deltaX, deltaY, beta
    double precision :: Re = 400, deltaT = 0.001

    deltaX = xlen/xGrid
    deltaY = ylen/yGrid

    !Create a staggered grid for u, v and p both for new and old value
    allocate(u(xGrid+1,yGrid+2))
    allocate(v(xGrid+2,yGrid+1))
    allocate(p(xGrid+2, yGrid+2))
    allocate(uNew(xGrid+1,yGrid+2))
    allocate(vNew(xGrid+2, yGrid+1))
    allocate(pNew(xGrid+2, yGrid+2))
    allocate(u2(xGrid+2, yGrid+2))
    allocate(v2(xGrid+2,yGrid+2))
    allocate(uv(xGrid+1,yGrid+1))
    allocate(rhsu(xGrid+1,yGrid+2))
    allocate(rhsv(xGrid+2,yGrid+1))
    allocate(psi(xGrid+1,yGrid+1))

    u = 1
    v = 1
    p = 1
    uNew = 0
    vNew = 0
    pNew = 0
    rhsu = 0
    rhsv = 0
    u2=0
    v2=0
    uv=0

    !Set the boudnary condition for variables that are on the boudnary
    uNew(1, :) = 0 !no slip at left wall
    uNew(xGrid+1, :) = 0 !no slip at right wall
    vNew(:, 1) = 0 !no slip at bottom wall
    vNew(:,yGrid+1) = 0 !Boundary condition at top wall v=0

    do WHILE(.not.(sum(abs(u - uNew))/((xGrid+1)*(yGrid+2)) < tolerance .and. & 
        sum(abs(v - vNew))/((xGrid+2)*(yGrid+1)) < tolerance))

        !Set boundary condition that are not coinciding with boundary condition
        uNew(1:xGrid+1,1) = -uNew(2:xGrid+2,2) ! no slip at bottom wall
        uNew(1:xGrid+1,yGrid+2) =2*uVelTop -uNew(1:xGrid+1, yGrid+1) !Boundary condition at the top wall
        vNew(1,1:yGrid+1) =  -vNew(2,1:yGrid+1) !no slip condition at left wall
        vNew(xGrid+2,1:yGrid+1) = -vNew(xGrid+1,1:yGrid+1) !no slip right wall

        u = uNew !update the u velocity
        v = vNew

        u2(2:xGrid+1,1:yGrid+2)=1.0/4*(u(2:xGrid+1,1:yGrid+2)+u(1:xGrid,1:yGrid+2))**2
        v2(1:xGrid+2,2:yGrid+1)=1.0/4*(v(1:xGrid+2,2:yGrid+1)+v(1:xGrid+2,1:yGrid))**2
        uv(1:xGrid+1,1:yGrid+1)=1.0/4*(u(1:xGrid+1,2:yGrid+2)+u(1:xGrid+1,1:yGrid+1))*&
                                (v(2:xGrid+2,1:yGrid+1)+v(1:xGrid+1,1:yGrid+1))
        
        rhsu(2:xGrid,2:yGrid+1)=u(2:xGrid,2:yGrid+1)-deltaT/deltaX*(u2(3:xGrid+1,2:yGrid+1)-u2(2:xGrid,2:yGrid+1))-&
                deltaT/deltaY*(uv(2:xGrid,2:yGrid+1)-uv(2:xGrid,1:yGrid))+deltaT/(Re*deltaX**2)*(u(1:xGrid-1,2:yGrid+1)-&
                2*u(2:xGrid,2:yGrid+1)+u(3:xGrid+1,2:yGrid+1))+deltaT/(Re*deltaY**2)*(u(2:xGrid,1:yGrid)-&
                2*u(2:xGrid,2:yGrid+1)+u(2:xGrid,3:yGrid+2))

        rhsu(1:xGrid+1,1) = -rhsu(2:xGrid+2,2) ! no slip at bottom wall
        rhsu(1:xGrid+1,yGrid+2) =2*uVelTop -rhsu(1:xGrid+1, yGrid+1) !Boundary condition at the top wall                        

        rhsv(2:xGrid+1,2:yGrid)=v(2:xGrid+1,2:yGrid)-deltaT/deltaY*(v2(2:xGrid+1,3:yGrid+1)-v2(2:xGrid+1,2:yGrid))-&
                deltaT/deltaX*(uv(2:xGrid+1,2:yGrid)-uv(1:xGrid,2:yGrid))+deltaT/(Re*deltaX**2)*(v(1:xGrid,2:yGrid)-&
                2*v(2:xGrid+1,2:yGrid)+v(3:xGrid+2,2:yGrid))+deltaT/(Re*deltaY**2)*(v(2:xGrid+1,1:yGrid-1)-&
                2*v(2:xGrid+1,2:yGrid)+v(2:xGrid+1,3:yGrid+1))

        rhsv(1,1:yGrid+1) =  -rhsv(2,1:yGrid+1) !no slip condition at left wall
        rhsv(xGrid+2,1:yGrid+1) = -rhsv(xGrid+1,1:yGrid+1) !no slip right wall

        beta=deltaX/deltaY
        p=1
                  
        do while(sum(abs(p - pNew))/((xGrid)*(yGrid)) > tolerance)
        !write(*,*) sum(abs(p - pNew))/((xGrid)*(yGrid))
        pNew(1,2:yGrid+1) = pNew(2,2:yGrid+1) !Presure normal is zero at bottom wall
        pNew(2:xGrid+1,1) = pNew(2:xGrid+1,2) !at left wall 
        pNew(xGrid+2,2:yGrid+1) = pNew(xGrid+1,2:yGrid+1) !at top wall
        pNew(2:xGrid+1,yGrid+2) = pNew(2:xGrid+1,yGrid+1) !at right wall
        p = pNew !update the pressure value

        pNew(2:xGrid+1,2:yGrid+1) = 1.0/(2*(1+beta**2))*(p(3:xGrid+2,2:yGrid+1)+p(1:xGrid,2:yGrid+1)+&
                            beta**2*(p(2:xGrid+1,3:yGrid+2)+p(2:xGrid+1,1:yGrid))-&
                            deltaX**2/deltaT*((rhsu(2:xGrid+1,2:yGrid+1)-rhsu(1:xGrid,2:yGrid+1))/deltaX+&
                            (rhsv(2:xGrid+1,2:yGrid+1)-rhsv(2:xGrid+1,1:yGrid))/deltaY))

        end do
    
        write(*,*) sum(abs(p - pNew))/((xGrid)*(yGrid))

        uNew(2:xGrid,2:yGrid+1)=-deltaT/deltaX*(pNew(3:xGrid+1,2:yGrid+1)-pNew(2:xGrid,2:yGrid+1))+rhsu(2:xGrid,2:yGrid+1)

        !find v at next step
        vNew(2:xGrid+1,2:yGrid)=-deltaT/deltaY*(pNew(2:xGrid+1,3:yGrid+1)-pNew(2:xGrid+1,2:xGrid))+rhsv(2:xGrid+1,2:yGrid)

        write(*,*) sum(abs(u - uNew))/((xGrid+1)*(yGrid+2)),sum(abs(v - vNew))/((xGrid+2)*(yGrid+1))
        
        if ( sum(abs(u - uNew))/((xGrid+1)*(yGrid+2)) > 10e6 ) then
            exit
        end if
        !exit 
    end do

    psi(1,:)=0

    do i=2,yGrid+1
        psi(i,:)=psi(i-1,:)+deltaX*v(i,:)
    end do

    OPEN(1, FILE = "allLid.m", STATUS='unknown',recl=2500) !contour of presssure
    
    write(1,*) "p=["
    do i=2,xGrid+1
        write(1,'(2x,128f16.12)') p(i,2:yGrid+1)
    end do
	write(1,*) "];"
	write(1,*)
	
	write(1,*) "u=["
    do i=2,xGrid+1
    	write(1,'(2x,128f16.12)') u(i,2:yGrid+1)
    end do
	write(1,*) "];"
	write(1,*)
	
    write(1,*) "v=["
    do i=2,xGrid+2
    	write(1,'(2x,128f16.12)') v(i,2:yGrid+1)
    end do
	write(1,*) "];"
	write(1,*)
	
	write(1,*) "psi=["
    do i=2,xGrid+2
    	write(1,'(2x,128f16.12)') psi(i,2:yGrid+1)
    end do
	write(1,*) "];"
	write(1,*)
	
    close(1)
    
    OPEN(2, FILE = "u.dat", STATUS='unknown') !contour of u velocity
    
    do i=1,xGrid+1
        do j=2,yGrid+1
            write(2,*) xLen/xGrid*(i-1), yLen/yGrid*(j-1)-0.5*deltaY, u(i,j)
        end do
        write(2,*) 
    end do

    close(2)

    OPEN(3, FILE = "v.dat", STATUS='unknown') !contour of v velocity
    
    do i=2,xGrid+1
        do j=1,yGrid+1
            write(3,*) xLen/xGrid*(i-1)-0.5*deltaX, yLen/yGrid*(j-1), v(i,j)
        end do
        write(3,*) 
    end do

    close(3)

    OPEN(4, FILE = "psi.dat", STATUS='unknown') !contour of v velocity
    
    do i=1,xGrid+1
        do j=1,yGrid+1
            write(4,*) xLen/xGrid*(i-1), yLen/yGrid*(j-1), psi(i,j)
        end do
        write(4,*) 
    end do

    close(4)

    OPEN(5, FILE = "p.dat", STATUS='unknown') !contour of p velocity
    
    do i=1,xGrid+1
        do j=2,yGrid+1
            write(5,*) xLen/xGrid*(i-1)-0.5*deltaX, yLen/yGrid*(j-1)-0.5*deltaY, p(i,j)
        end do
        write(5,*) 
    end do

    close(5)
    
    OPEN(6, FILE = "Vel.dat", STATUS='unknown', recl=200) !contour of u velocity
    
    do i=1,xGrid+1
        do j=1,yGrid+1
            write(6,*) xLen/xGrid*(i-1), yLen/yGrid*(j-1), u(i,j), v(i,j)
        end do
        write(6,*) 
    end do

    close(6)

end program MAC

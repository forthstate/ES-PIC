C     Utilities file containing some expensive functions: Gather,Scatter

C     loop to interpolate field-data from surrounding nodes to a given position
      subroutine gatherEF(lc,efx,efy,efz,nx,ny,nz,fx,fy,fz)
      integer nx, ny, nz
      real*8 efx(0:nx-1, 0:ny-1, 0:nz-1), efy(0:nx-1, 0:ny-1, 0:nz-1)
      real*8 efz(0:nx-1, 0:ny-1, 0:nz-1), lc(0:2)
      real*8 fx, fy, fz
Cf2py intent(in) efx, efy, efz
Cf2py intent(out) fx, fy, fz
      integer i, j, k 
      
      
C      lc(0) = (x(0) - x0(0))/dx
      i = int(lc(0)); di = lc(0)-i
C      lc(1) = (x(1) - x0(1))/dy
      j = int(lc(1)); dj = lc(1)-j
C      lc(2) = (x(2) - x0(2))/dz
      k = int(lc(2)); dk = lc(2)-k
            
      fx = efx(i,j,k)*(1-di)*(1-dj)*(1-dk) +
     &     efx(i+1,j,k)*(di)*(1-dj)*(1-dk) +
     &     efx(i,j+1,k)*(1-di)*(dj)*(1-dk) +
     &     efx(i,j,k+1)*(1-di)*(1-dj)*(dk) +
     &     efx(i+1,j+1,k)*(di)*(dj)*(1-dk) +
     &     efx(i,j+1,k+1)*(1-di)*(dj)*(dk) +
     &     efx(i+1,j,k+1)*(di)*(1-dj)*(dk) +
     &     efx(i+1,j+1,k+1)*(di)*(dj)*(dk)
     
      fy = efy(i,j,k)*(1-di)*(1-dj)*(1-dk) +
     &     efy(i+1,j,k)*(di)*(1-dj)*(1-dk) +
     &     efy(i,j+1,k)*(1-di)*(dj)*(1-dk) +
     &     efy(i,j,k+1)*(1-di)*(1-dj)*(dk) +
     &     efy(i+1,j+1,k)*(di)*(dj)*(1-dk) +
     &     efy(i,j+1,k+1)*(1-di)*(dj)*(dk) +
     &     efy(i+1,j,k+1)*(di)*(1-dj)*(dk) +
     &     efy(i+1,j+1,k+1)*(di)*(dj)*(dk)
     
      fz = efz(i,j,k)*(1-di)*(1-dj)*(1-dk) +
     &     efz(i+1,j,k)*(di)*(1-dj)*(1-dk) +
     &     efz(i,j+1,k)*(1-di)*(dj)*(1-dk) +
     &     efz(i,j,k+1)*(1-di)*(1-dj)*(dk) +
     &     efz(i+1,j+1,k)*(di)*(dj)*(1-dk) +
     &     efz(i,j+1,k+1)*(1-di)*(dj)*(dk) +
     &     efz(i+1,j,k+1)*(di)*(1-dj)*(dk) +
     &     efz(i+1,j+1,k+1)*(di)*(dj)*(dk)
      
      return
      end
      
C     loop to extrapolate particle-data to surrounding nodes
      subroutine scatterMPW(lc,den,mpw,nx,ny,nz)
      integer nx,ny,nz
      real*8 den(0:nx-1, 0:ny-1, 0:nz-1)
      real*8 mpw, lc(0:2)
Cf2py intent(in, out) den
      integer i, j, k 
      
C     lc(0) = (x(0) - x0(0))/dx
      i = int(lc(0)); di = lc(0)-i
C     lc(1) = (x(1) - x0(1))/dy
      j = int(lc(1)); dj = lc(1)-j
C     lc(2) = (x(2) - x0(2))/dz
      k = int(lc(2)); dk = lc(2)-k
      
      den(i,j,k)=den(i,j,k)+mpw*(1-di)*(1-dj)*(1-dk)
      den(i+1,j,k)=den(i+1,j,k)+mpw*(di)*(1-dj)*(1-dk)
      den(i,j+1,k)=den(i,j+1,k)+mpw*(1-di)*(dj)*(1-dk)
      den(i,j,k+1)=den(i,j,k+1)+mpw*(1-di)*(1-dj)*(dk)
      den(i+1,j+1,k)=den(i+1,j+1,k)+mpw*(di)*(dj)*(1-dk)
      den(i+1,j,k+1)=den(i+1,j,k+1)+mpw*(di)*(1-dj)*(dk)
      den(i,j+1,k+1)=den(i,j+1,k+1)+mpw*(1-di)*(dj)*(dk)
      den(i+1,j+1,k+1)=den(i+1,j+1,k+1)+mpw*(di)*(dj)*(dk)
      
      return
      end

C     loop to advance particles
      subroutine advancePAR(lc,efx,efy,efz,nx,ny,nz,v,p,r1,r2,dt,dm)
      integer nx, ny, nz
      real*8 efx(0:nx-1, 0:ny-1, 0:nz-1), efy(0:nx-1, 0:ny-1, 0:nz-1)
      real*8 efz(0:nx-1, 0:ny-1, 0:nz-1), lc(0:2), v(0:2), p(0:2)
      real*8 fx, fy, fz, r1(0:2), r2(0:2), dt, dm
Cf2py intent(in) efx, efy, efz
Cf2py intent(in, out) v, p
      integer i, j, k
      
      i = int(lc(0)); di = lc(0)-i
      j = int(lc(1)); dj = lc(1)-j
      k = int(lc(2)); dk = lc(2)-k
            
      fx = efx(i,j,k)*(1-di)*(1-dj)*(1-dk) +
     &     efx(i+1,j,k)*(di)*(1-dj)*(1-dk) +
     &     efx(i,j+1,k)*(1-di)*(dj)*(1-dk) +
     &     efx(i,j,k+1)*(1-di)*(1-dj)*(dk) +
     &     efx(i+1,j+1,k)*(di)*(dj)*(1-dk) +
     &     efx(i,j+1,k+1)*(1-di)*(dj)*(dk) +
     &     efx(i+1,j,k+1)*(di)*(1-dj)*(dk) +
     &     efx(i+1,j+1,k+1)*(di)*(dj)*(dk)
     
      fy = efy(i,j,k)*(1-di)*(1-dj)*(1-dk) +
     &     efy(i+1,j,k)*(di)*(1-dj)*(1-dk) +
     &     efy(i,j+1,k)*(1-di)*(dj)*(1-dk) +
     &     efy(i,j,k+1)*(1-di)*(1-dj)*(dk) +
     &     efy(i+1,j+1,k)*(di)*(dj)*(1-dk) +
     &     efy(i,j+1,k+1)*(1-di)*(dj)*(dk) +
     &     efy(i+1,j,k+1)*(di)*(1-dj)*(dk) +
     &     efy(i+1,j+1,k+1)*(di)*(dj)*(dk)
     
      fz = efz(i,j,k)*(1-di)*(1-dj)*(1-dk) +
     &     efz(i+1,j,k)*(di)*(1-dj)*(1-dk) +
     &     efz(i,j+1,k)*(1-di)*(dj)*(1-dk) +
     &     efz(i,j,k+1)*(1-di)*(1-dj)*(dk) +
     &     efz(i+1,j+1,k)*(di)*(dj)*(1-dk) +
     &     efz(i,j+1,k+1)*(1-di)*(dj)*(dk) +
     &     efz(i+1,j,k+1)*(di)*(1-dj)*(dk) +
     &     efz(i+1,j+1,k+1)*(di)*(dj)*(dk)
      
C     update velocity
      v(0) = v(0) + fx*dm
      v(1) = v(1) + fy*dm
      v(2) = v(2) + fz*dm
C     update position
      p(0) = p(0) + v(0)*dt
      p(1) = p(1) + v(1)*dt
      p(2) = p(2) + v(2)*dt
      
C     REFLECT Particles at the Walls
      do l=0,2
        if(p(l).lt.r1(l))then
          p(l) = 2*r1(l) - p(l)
          v(l) = -1*v(l)
        else if(p(l).ge.r2(l))then
          p(l) = 2*r2(l) - p(l)
          v(l) = -1*v(l)
        end if
      end do
      
      return
      end    
      
      
            

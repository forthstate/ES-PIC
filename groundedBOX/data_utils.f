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
      real*8 di, dj, dk
      
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
      
      return
      end
      
C     loop to extrapolate particle-data to surrounding nodes
      subroutine scatterMPW(lc,den,mpw,nx,ny,nz)
      integer nx,ny,nz
      real*8 den(0:nx-1, 0:ny-1, 0:nz-1), lc(0:2)
      real*8 mpw
Cf2py intent(in, out) den
      integer i, j, k 
      real*8 di, dj, dk
      
      i = int(lc(0)); di = lc(0)-i
      j = int(lc(1)); dj = lc(1)-j
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







      

C     Utilities file containing some expensive functions: Potential,LeastSquareError,ElectricField

C     loop to calcalute the potential at nodes
      subroutine potential_loop(rho,phi,nx,ny,nz,idx2,idy2,idz2,EP0,w)
      integer nx, ny, nz
      real*8 phi(0:nx-1, 0:ny-1, 0:nz-1), rho(0:nx-1, 0:ny-1, 0:nz-1)
      real*8 idx2, idy2, idz2, EP0, w
Cf2py intent(in, out) phi
Cf2py intent(in) rho

      do i=1,nx-2
        do j=1,ny-2
          do k=1,nz-2    
            phi(i,j,k) = (1-w)*phi(i,j,k) + w*(rho(i,j,k)/EP0 + 
     &          idx2*(phi(i-1,j,k) + phi(i+1,j,k)) +
     &          idy2*(phi(i,j-1,k) + phi(i,j+1,k)) + 
     &          idz2*(phi(i,j,k-1) + phi(i,j,k+1))) /
     &          (2*idx2 + 2*idy2 + 2*idz2)
          end do
        end do
      end do
      return
      end

C     loop get least squared array for convergence check      
      subroutine convergence_loop(rho,phi,nx,ny,nz,idx2,idy2,idz2,EP0,s)
      integer nx, ny, nz
      real*8 phi(0:nx-1, 0:ny-1, 0:nz-1), rho(0:nx-1, 0:ny-1, 0:nz-1)
      real*8 idx2, idy2, idz2, EP0, s
Cf2py intent(in) rho, phi      
Cf2py intent(out) s
      real*8 R

      s = 0
      do i=1,nx-2
        do j=1,ny-2
          do k=1,nz-2    
            R = -phi(i,j,k)*(2*idx2 + 2*idy2 + 2*idz2) + 
     &          rho(i,j,k)/EP0 +             
     &          idx2*(phi(i-1,j,k) + phi(i+1,j,k)) +
     &          idy2*(phi(i,j-1,k) + phi(i,j+1,k)) + 
     &          idz2*(phi(i,j,k-1) + phi(i,j,k+1))
     
          s = s + R*R
          end do
        end do
      end do
      return
      end
      
C     loop to calculate electric field      
      subroutine electric_field(phi,efx,efy,efz,nx,ny,nz,dx,dy,dz)
      integer nx, ny, nz
      real*8 phi(0:nx-1, 0:ny-1, 0:nz-1), efx(0:nx-1, 0:ny-1, 0:nz-1)
      real*8 efy(0:nx-1, 0:ny-1, 0:nz-1), efz(0:nx-1, 0:ny-1, 0:nz-1)
      real*8 dx, dy, dz
Cf2py intent(in, out) efx, efy, efz
Cf2py intent(in) phi

      do i=0,nx-1
        do j=0,ny-1
          do k=0,nz-1
          
C           x-component
            if(i.eq.0)then
              efx(i,j,k) = -(-3*phi(i,j,k) + 4*phi(i+1,j,k) +
     &                     phi(i+2,j,k))/(2*dx)
            else if(i.eq.nx-1)then
              efx(i,j,k)= -(phi(i-2,j,k) - 4*phi(i-1,j,k) + 
     &                    3*phi(i,j,k))/(2*dx)
            else
              efx(i,j,k) = -(phi(i+1,j,k) - phi(i-1,j,k))/(2*dx)
            end if              
C           y-component
            if(j.eq.0)then
              efy(i,j,k) = -(-3*phi(i,j,k) + 4*phi(i,j+1,k) -
     &                     phi(i,j+2,k))/(2*dy)              
            else if(j.eq.ny-1)then
              efy(i,j,k)= -(phi(i,j-2,k) - 4*phi(i,j-1,k) +
     &                    3*phi(i,j,k))/(2*dy)              
            else
              efy(i,j,k) = -(phi(i,j+1,k) - phi(i,j-1,k))/(2*dy)
            end if
C           z-component
            if(k.eq.0)then
              efz(i,j,k) = -(-3*phi(i,j,k) + 4*phi(i,j,k+1) -
     &                     phi(i,j,k+2))/(2*dz)              
            else if(k.eq.nz-1)then
              efz(i,j,k)= -(phi(i,j,k-2) - 4*phi(i,j,k-1) +
     &                    3*phi(i,j,k))/(2*dz)              
            else
              efz(i,j,k) = -(phi(i,j,k+1) - phi(i,j,k-1))/(2*dz)
            end if
                                                     
          end do
        end do
      end do
      return
      end
            
      
      
      
            


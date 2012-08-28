 subroutine box_qp_f(Q,u,b,rho,MAXIT,tol,qq,grad)
      implicit none
      integer :: qq, MAXIT, outer, j
!     double precision Q(qq,qq), u(qq), b(qq), grad(qq)
!     double precision tol, tt, bb, uold, rho, objcur, del, objold, dlx 
      real(kind=8), dimension (qq,qq) :: Q
      real(kind=8), dimension (qq) :: u, b, grad
      real(kind=8) :: tol, tt, bb, uold, rho, objcur, del, objold, dlx 
      grad = MATMUL(Q,u+b)
      grad=2.0*grad
!      write(*,*) grad(1:2)
      objcur=DOT_PRODUCT(grad, b+u)
      objold=objcur
      outer=0
      dlx=0.0
10231 continue
      outer=outer+1    
10240 do 10241 j=1,qq                                                        
      uold=u(j)   
!     write(*,*) uold  
      bb=Q(j,j)*u(j)                                                          
      bb=grad(j) - (2.0*bb) 
      tt=abs(bb)/(2.0*Q(j,j)) 
      if(tt.lt.rho) then 
        u(j)=sign(tt,-bb) 
      else 
         u(j)=sign(rho,-bb)  
      end if
      if(u(j).eq.uold) goto 10241                                              
      del=u(j)- uold       
      grad=grad + 2.0*del*Q(:,j)  
!      write(*,*) grad(1:2)                                                     
10241 continue                                                              
      objcur=DOT_PRODUCT(grad, b+u)
      dlx = abs(objcur - objold)/(abs(objold)+0.000001) 
      objold=objcur
!     write(*,*) (dlx.lt.tol).or.(outer.gt.(MAXIT-1)), objcur, dlx, (dlx.lt.tol)
      if ((dlx.lt.tol).or.(outer.gt.(MAXIT-1))) goto 10232    
      goto 10231                                                            
10232 continue  
10233 continue                                                            
      return
      end  
                                                






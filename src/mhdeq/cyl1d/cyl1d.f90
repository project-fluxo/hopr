!============================================================================================================ xX =================
!        _____     _____    _______________    _______________   _______________            .xXXXXXXXx.       X
!       /    /)   /    /)  /    _____     /)  /    _____     /) /    _____     /)        .XXXXXXXXXXXXXXx  .XXXXx
!      /    //   /    //  /    /)___/    //  /    /)___/    // /    /)___/    //       .XXXXXXXXXXXXXXXXXXXXXXXXXx
!     /    //___/    //  /    //   /    //  /    //___/    // /    //___/    //      .XXXXXXXXXXXXXXXXXXXXXXX`
!    /    _____     //  /    //   /    //  /    __________// /    __      __//      .XX``XXXXXXXXXXXXXXXXX`
!   /    /)___/    //  /    //   /    //  /    /)_________) /    /)_|    |__)      XX`   `XXXXX`     .X`
!  /    //   /    //  /    //___/    //  /    //           /    //  |    |_       XX      XXX`      .`
! /____//   /____//  /______________//  /____//           /____//   |_____/)    ,X`      XXX`
! )____)    )____)   )______________)   )____)            )____)    )_____)   ,xX`     .XX`
!                                                                           xxX`      XXx
! Copyright (C) 2017  Florian Hindenlang <hindenlang@gmail.com>
! This file is part of HOPR, a software for the generation of high-order meshes.
!
! HOPR is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! HOPR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with HOPR. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "hopr.h"
MODULE MOD_Cyl1d
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitCyl1d 
  MODULE PROCEDURE InitCyl1d 
END INTERFACE

! allow different dimensions of input/output arrays
!INTERFACE MapToCyl1d 
!  MODULE PROCEDURE MapToCyl1d 
!END INTERFACE

INTERFACE FinalizeCyl1d 
  MODULE PROCEDURE FinalizeCyl1d 
END INTERFACE

PUBLIC::InitCyl1d
PUBLIC::MapToCyl1d
PUBLIC::FinalizeCyl1d
!===================================================================================================================================

CONTAINS

SUBROUTINE InitCyl1d 
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals, ONLY:UNIT_StdOut
USE MOD_ReadInTools, ONLY: GETINT,GETREAL,GETREALARRAY
USE MOD_Cyl1d_Vars
!USE MOD_CCint,ONLY:InitCCint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'  INIT CYL1D INPUT ...'

Bz0=GETREAL('cyl1d_Bz0',"1.")
alpha0=GETREAL('cyl1d_alpha0') !TOK: 0.625, PP: 4.
eta_param=GETREALARRAY('cyl1d_eta_param',3)  !TOK: (24,4,1.5), PP:(20,10,1)

nps=100 !number of equidistant points used for the spline
ALLOCATE(rSpl(0:nps),ohmBt_spl(4,0:nps),ohmBz_Spl(4,0:nps))

!initialize Cyl1d parameters
CALL ohm_equil_1d(10)

WRITE(UNIT_stdOut,'(A)')'  ... DONE'
END SUBROUTINE InitCyl1d

!=================================================================================================================================
!> 1D equilibrium code, by D. Bonfiglio (used in Specyl and Pixie3d for initialization) 
!> using profile parameter alpha0=(E_0/eta_0), with radially varying eta=eta_0*(1+Ar^B)^C
!> iteration of Bz and Bt via integral relation from Ohms law and force balance
!=================================================================================================================================
SUBROUTINE ohm_equil_1d(n_int)
! MODULES
USE MOD_Globals
USE MOD_Cyl1d_vars ,ONLY: nps,alpha0,Bz0,rSpl,ohmBt_Spl,ohmBz_Spl
USE SPLINE1_MOD    ,ONLY: SPLINE1_FIT,SPLINE1_EVAL
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
  INTEGER, INTENT(IN) :: n_int   !< number of points for the integrati
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                      :: np           !number of total grid points
  INTEGER                      :: ix,it
  INTEGER                      :: ioUnit=999
  INTEGER                      :: iGuess=1
  REAL                         :: splout(3)
  REAL                         :: bzh,bth,dr,err
  REAL,DIMENSION(0:nps*n_int)  :: ohmBt,ohmBz,vr_over_eta0  !result of iteration
  REAL,DIMENSION(0:nps*n_int)  :: rg,ohmBt_old,ohmBz_old,dummy
!  REAL,DIMENSION(0:nps*n_int-1):: rh,etah 
  REAL,DIMENSION(0:nps*n_int)  :: eta,eta_b2 
  LOGICAL                      :: converged
!=================================================================================================================================
WRITE(UNIT_stdOut,'(A)')' 1D OHMIC equilibrium solver, constant pressure'
  np=nps*n_int !number of total grid points

  !initialize guess for iterative procedure
  ohmBz(0:np)=Bz0
  !Initial guess for B_theta
  dr=1./REAL(np)
  rg(0)=0.
  ohmBt(0)=0.
  eta(0)= eta_rprofile(rg(0))
  DO ix=1,np
    rg(ix)    = REAL(ix)*dr          !radial grid [0,1]
    ohmBt(ix) = 0.5*rg(ix)*alpha0
!    rh(ix-1)  = 0.5*(rg(ix)-rg(ix-1)) !half grid
!    etah(ix-1)= eta_rprofile(rh(ix-1))
    eta(ix)= eta_rprofile(rg(ix))
  END DO

  rSpl(0)=0.
  DO ix=1,nps
    rSpl(ix)    = REAL(ix)/REAL(nps) !radial grid for spline [0,1]
  END DO


  converged=.FALSE.
  !Iteration
  DO it=1,100

    ohmBt_old(0:np) = ohmBt(0:np)
    ohmBz_old(0:np) = ohmBz(0:np)

    eta_b2(:)=eta(:)*(ohmbt(:)**2+ohmbz(:)**2)
    !Perform integral of B_theta
    dummy(0) = 0.0d0
    ohmBt(0) = 0.0d0
    DO ix=1,np

!      bzh = 0.5*(ohmBz(ix)+ohmBz(ix-1))
!      bth = 0.5*(ohmBt(ix)+ohmBt(ix-1))
!
!      dummy(ix) = dummy(ix-1) + & 
!                   bzh**2/(etah(ix-1)*(bzh**2+bth**2))*rh(ix-1)*dr
      dummy(ix) = dummy(ix-1) + & 
                   dr*0.5*( ohmbz(ix  )**2*rg(ix  )/eta_b2(ix  ) &
                           +ohmbz(ix-1)**2*rg(ix-1)/eta_b2(ix-1) )
      !WRITE(*,*) rh(ix-1),etah(ix-1))

      ohmBt(ix) = alpha0*dummy(ix)/rg(ix)

      !WRITE(*,*) ' ix, ohmBt(ix): ', ohmBt(ix)

    END DO !ix=1,np 

    eta_b2(:)=eta(:)*(ohmbt(:)**2+ohmbz(:)**2)
    !Perform integral of Bz
    dummy(0) = 0.0d0
    DO ix=1,np
!      bzh = 0.5*(ohmBz(ix)+ohmBz(ix-1))
!      bth = 0.5*(ohmBt(ix)+ohmBt(ix-1))
!
!      dummy(ix) = dummy(ix-1) + bth/(etah(ix-1)*(bzh**2+bth**2))*dr 

      dummy(ix) = dummy(ix-1) + dr*0.5*( ohmbt(ix  )/eta_b2(ix  ) &
                                        +ohmbt(ix-1)/eta_b2(ix-1) )
           
      !WRITE(*,*) rh(ix-1),etah(ix-1)

      ohmBz(ix) = Bz0*EXP(-alpha0*dummy(ix))

      !WRITE(*,*)  ' ix, ohmBz(ix): ', ohmBz(ix)

    END DO !ix=1,np 

    !interpolate to spline with nps points and reevaluate at np points
    ohmBt_Spl=0.
    ohmBz_spl=0.
    ohmBt_Spl(1,0:nps)=ohmBt(0:np:n_int)
    ohmBz_spl(1,0:nps)=ohmBz(0:np:n_int)
    CALL SPLINE1_FIT(nps+1,rSpl,ohmBt_Spl(:,:), K_BC1=4, K_BCN=0)
    CALL SPLINE1_FIT(nps+1,rSpl,ohmBz_Spl(:,:), K_BC1=3, K_BCN=0) !bz'=0 at r=0
    DO ix=1,np
      CALL SPLINE1_EVAL((/1,0,0/), nps+1,rg(ix),rSpl,ohmBt_Spl(:,:),iGuess,splout) 
      ohmBt(ix)=splout(1)
      CALL SPLINE1_EVAL((/1,0,0/), nps+1,rg(ix),rSpl,ohmBz_Spl(:,:),iGuess,splout) 
      ohmBz(ix)=splout(1)
    END DO
    
    !Check convergence
    err = SUM((ohmBt(:)-ohmBt_old(:))**2+ (ohmBz(:)-ohmBz_old(:))**2 )
    err = sqrt(0.5*err*dr)

    !WRITE(*,*) 'Equilibrium iter =',it,' ; Error=',err

    converged= (err < 1.0e-6*dr**2)
    
    IF(converged) EXIT

  END DO !it=1,100 

  vr_over_eta0(:) = -alpha0*ohmBt(:)/(ohmBz(:)**2.+ohmBt(:)**2.)


  OPEN(UNIT   = ioUnit       ,&
       FILE   = TRIM(Projectname)//'_cyl1d.csv'  ,&
       STATUS = 'Unknown'    ,&
       ACCESS = 'SEQUENTIAL' ) 
    WRITE(ioUnit,'(A)')' "r", "btheta", "bz", "vr_over_eta0", "eta"'
    DO ix=0,np
      WRITE(ioUnit,'(*(E23.15,:,","))') rg(ix),ohmBt(ix),ohmBz(ix),vr_over_eta0(ix),eta_rprofile(rg(ix))
    END DO
  CLOSE(ioUnit)

  WRITE(*,*)
  IF(converged)THEN
    WRITE(*,'(a,i3,a,e10.2)') ' Equilibrium converged in ',it &
                              ,' iterations with error =',err
     
  END IF


  IF(.NOT.converged) STOP '1d equilibrium solver not converged!!'

WRITE(UNIT_stdOut,'(A)')' ...DONE.'
END SUBROUTINE ohm_equil_1d


!=================================================================================================================================
!> compute the radially dependant resistivity, scaling the constant resistivity  eta/eta_0=(1+A*r^B)^C
!> needs initialization parameters eta_param=(A,B,C)
!=================================================================================================================================
PURE FUNCTION eta_rprofile(rpos)
! MODULES
USE MOD_cyl1d_vars ,ONLY: eta_param !profile parameters: A,B,C
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: rpos  !< radial position [0,1]
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: eta_rprofile !< eta at radial position
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VaLIABLES
!==================================================================================================================================
eta_rprofile=(1.0+eta_param(1)*(rpos**eta_param(2)))**eta_param(3)
END FUNCTION eta_rprofile

SUBROUTINE MapToCyl1d(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
!===================================================================================================================================
! Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from a tokamak soloviev equilibrium. 
! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
! for a fixed theta, uses a newton method to find the exact location  in r of psi(x(r),y(r))-Psi0=0, psi0=psi(psinorm=r_p^2),
! using an approximate map(R,Z)<->(rho,theta) of the soloviev equilibrium. 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars ,ONLY: nVarMHDEQ
USE MOD_MHDEQ_Vars ,ONLY: nRhoCoefs,RhoCoefs
USE MOD_MHDEQ_Tools,ONLY: Eval1DPoly
USE MOD_Cyl1d_Vars ,ONLY: nps,rSpl,ohmBt_Spl,ohmBz_Spl,alpha0
USE SPLINE1_MOD    ,ONLY: SPLINE1_EVAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nTotal         ! total number of points
REAL, INTENT(IN)   :: x_in(3,nTotal) ! input coordinates represent a cylinder: 
INTEGER, INTENT(IN):: InputCoordSys  ! =0: x_in(1:3) are (x,y,z) coordinates in a cylinder of size r=[0;1], z=[0;1]
                                     ! =1: x_in(1:3) are (r,z,phi) coordinates r= [0;1], z= [0;1], phi=[0;1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: x_out(3,nTotal) ! mapped x,y,z coordinates with vmec data
REAL,INTENT(OUT)   :: MHDEQdata(nVarMHDEQ,nTotal) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: r_p    ! raduis in cylindrical coordinate system
REAL    :: theta  ! some poloidal angle [0,2pi] ,THIS IS NOT =atan((Z-Zaxis)/(R-Raxis))
INTEGER :: iNode
INTEGER :: percent
INTEGER :: iGuess=1
REAL    :: splout(3)
REAL    :: Density
REAL    :: Bt,Bz,vr_over_eta0,Bcart(3) 
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,I8,A,A,A)')'  MAP ', nTotal,' NODES TO CYLINDER EQUILIBRIUM'
percent=0
DO iNode=1,nTotal
  ! output of progress in %
  IF((nTotal.GT.10000).AND.(MOD(iNode,(nTotal/100)).EQ.0)) THEN
    percent=percent+1
    WRITE(0,'(I4,A23,A1)',ADVANCE='NO')percent, ' % of nodes evaluated...',ACHAR(13)
  END IF
  SELECT CASE(InputCoordSys)
  CASE(0)!x_in(1:3) = x,y,z of cylinder with r<1 and z=[0;1]
    r_p   = SQRT(x_in(1,iNode)**2+x_in(2,iNode)**2) 
    Theta = ATAN2(x_in(2,iNode),x_in(1,iNode))
    x_out(1:3,iNode)= x_in(1:3,iNode)
  CASE(1) !x_in(1:3) = (r,theta,z) with r= [0;1], theta=[0;1] 
    r_p =  x_in(1,iNode) !=r
    Theta = 2.*Pi*x_in(3,iNode) !=2*pi*phi
    x_out(1,iNode)= r_p*COS(theta)
    x_out(2,iNode)= r_p*SIN(theta)
    x_out(3,iNode)= x_in(3,iNode)
  END SELECT 
  IF(r_p.GT.1.0+1.0e-08) THEN
    WRITE(*,*) 'radius must be <=1',r_p
    STOP
  END IF

  Density=Eval1DPoly(nRhoCoefs,RhoCoefs,r_p) 

  CALL SPLINE1_EVAL((/1,0,0/),nps+1,r_p,rSpl,ohmBt_spl(:,:),iGuess,splout)
  Bt=splout(1)
  CALL SPLINE1_EVAL((/1,0,0/),nps+1,r_p,rSpl,ohmBz_spl(:,:),iGuess,splout)
  Bz=splout(1)
  vr_over_eta0 = -alpha0*Bt/(Bz**2.+Bt**2.)

  Bcart(1)=-Bt*SIN(theta)
  Bcart(2)= Bt*COS(theta)
  Bcart(3)= Bz

  MHDEQdata(:,iNode)=0.
  MHDEQdata(  1,iNode)=Density
  MHDEQdata(  2,iNode)=0.        !pressure not yet included
  MHDEQdata( 3:5,iNode)=Bcart(:)
  MHDEQdata(   6,iNode)=r_p !poloidal flux (dummy, normalized)
  MHDEQdata(   7,iNode)=r_p !toroidal flux (dummy, normalized)
  MHDEQdata(8:10,iNode)=0.
  MHDEQdata(11,iNode)=vr_over_eta0*COS(theta)  ! vx/eta_0 from radial velocity vr=vr_over_eta0*eta_0
  MHDEQdata(12,iNode)=vr_over_eta0*SIN(theta)  ! vy/eta_0
  MHDEQdata(13,iNode)=0.                       ! vz/eta_0
END DO !iNode


WRITE(UNIT_stdOut,'(A)')'  ...DONE.                             '


END SUBROUTINE MapToCyl1d 

!===================================================================================================================================
!> Finalize Sololviev module
!!
!===================================================================================================================================
SUBROUTINE FinalizeCyl1d
! MODULES
USE MOD_CCint,ONLY:FinalizeCCint
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL FinalizeCCint()


END SUBROUTINE FinalizeCyl1d


END MODULE MOD_Cyl1d

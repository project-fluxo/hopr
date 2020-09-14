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
MODULE MOD_GVEC
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

INTERFACE InitGVEC 
  MODULE PROCEDURE InitGVEC 
END INTERFACE

! allow different dimensions of input/output arrays
!INTERFACE MapToGVEC 
!  MODULE PROCEDURE MapToGVEC 
!END INTERFACE

INTERFACE FinalizeGVEC 
  MODULE PROCEDURE FinalizeGVEC 
END INTERFACE

PUBLIC::InitGVEC
PUBLIC::MapToGVEC
PUBLIC::FinalizeGVEC
!===================================================================================================================================

CONTAINS
SUBROUTINE InitGVEC 
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,abort
USE MOD_ReadInTools
USE MOD_GVEC_Vars
#ifdef PP_GVEC
USE MODgvec_gvec_to_hopr, ONLY: Init_gvec_to_hopr
#endif /*PP_GVEC*/
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
INTEGER                 :: GVEC_factorSFL
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'  INIT GVEC INPUT ...'

GVECdataFile=GETSTR("GVECwoutfile")

GVEC_SFLcoord=GETINT("GVEC_SFLcoord","0")
GVEC_factorSFL=GETINT("GVEC_factorSFL","4")

#ifdef PP_GVEC
CALL Init_gvec_to_hopr(GVECdataFile,SFLcoord_in=GVEC_SFLcoord,factorSFL_in=GVEC_factorSFL)
#else
STOP 'TRYING TO USE GVEC INTERFACE, BUT HOPR IS NOT LINKED TO GVEC!!!'
#endif /*PP_GVEC*/

WRITE(UNIT_stdOut,'(A)')'  ... DONE'
END SUBROUTINE InitGVEC

SUBROUTINE MapToGVEC(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
!===================================================================================================================================
! Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from GVEC data. 
! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars,    ONLY: nVarMHDEQ
USE MOD_MHDEQ_Vars,    ONLY: nRhoCoefs,RhoFluxVar,RhoCoefs
USE MOD_MHDEQ_Tools,   ONLY: Eval1DPoly
USE MOD_GVEC_Vars,     ONLY: wp,twoPi,mu0,GVECdatafile
#ifdef PP_GVEC
USE MODgvec_gvec_to_hopr , ONLY: gvec_to_hopr 
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nTotal          !! total number of points
REAL(wp),INTENT(IN):: x_in(3,nTotal)  !! input coordinates represent a cylinder: 
INTEGER, INTENT(IN):: InputCoordSys   !!  0: x_in(1:3) are (x,y,z) coordinates in a cylinder of size r=[0;1], z=[0;1]
                                      !! 1: x_in(1:3) are (r,z,phi) coordinates r= [0;1], z= [0;1], phi=[0;1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT):: x_out(3,nTotal) !! mapped x,y,z coordinates with GVEC data
REAL(wp),INTENT(OUT):: MHDEQdata(nVarMHDEQ,nTotal) !! vector of equilibrium variables, see definition in mhdeq_vars.f90
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode,minnode,maxnode,istep,nsteps,nnodes
REAL(wp) :: r_p, theta,zeta,phi_int,chi_int,phi_edge_axis(2),chi_edge_axis(2),Density,phinorm,chinorm
REAL(wp) :: xin2(3,nTotal)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,I8,A,A,A)')'  MAP ', nTotal,' NODES TO GVEC DATA FROM ',TRIM(GVECdataFile),' ...'
nnodes=nTotal/100
nsteps=nTotal/nNodes+1
CALL ProgressBar(0,nTotal)
DO istep=0,nsteps-1
  minNode=istep*nnodes+1
  IF(minNode.GT.nTotal) EXIT
  maxNode=MIN((istep+1)*nnodes,nTotal)
  DO iNode=minNode,maxNode
    SELECT CASE(InputCoordSys)
    CASE(0)!x_in(1:3) = x,y,z of cylinder with r<1 and z=[0;1]
      r_p   = SQRT(x_in(1,iNode)**2+x_in(2,iNode)**2) 
      theta = ATAN2(x_in(2,iNode),x_in(1,iNode))
      zeta  = -twoPi*x_in(3,iNode) 
    CASE(1) !x_in(1:3) = (r,z,phi) with r= [0;1], z= [0;1], phi=[0;1] 
      r_p =  x_in(1,iNode) ! =r
      theta = twoPi*x_in(3,iNode) ! =2*pi*phi
      zeta  = twoPi*x_in(2,iNode) ! =2*pi*z
    END SELECT 
  
    xin2(:,iNode)=(/r_p,theta,-zeta/)
  END DO

#ifdef PP_GVEC
  CALL gvec_to_hopr(maxnode-minnode+1,xin2(:,minnode:maxnode), &
                                  x_out(:,minnode:maxnode),&
                           MHDEQdata(2:10,minnode:maxnode),phi_edge_axis,chi_edge_axis)
#else
  x_out=0. 
  STOP 'gvec_to_hopr cannot be called, HOPR not linked to GVEC.'
#endif



  DO iNode=minNode,maxNode
    phi_int=MHDEQdata(7,iNode)
    chi_int=MHDEQdata(6,iNode)
    ! output of progress in %
    !NOTE: Phi is the toroidal flux, called PHI in GVEC
    Phinorm=(phi_int-phi_edge_axis(1))/(phi_edge_axis(2)-phi_edge_axis(1))
    Chinorm=(chi_int-chi_edge_axis(1))/(chi_edge_axis(2)-chi_edge_axis(1))
    
    SELECT CASE (RhoFluxVar)
    CASE(0) !use normalized toroidal flux 
      Density=Eval1DPoly(nRhoCoefs,RhoCoefs,Phinorm) 
    CASE(1) !use normalized poloidal flux 
      Density=Eval1DPoly(nRhoCoefs,RhoCoefs,chinorm) 
    CASE(2) !use sqrt of normalized toroidal flux 
      Density=Eval1DPoly(nRhoCoefs,RhoCoefs,sqrt(Phinorm))
    CASE(3) !use sqrt of normalized poloidal flux 
      Density=Eval1DPoly(nRhoCoefs,RhoCoefs,sqrt(chinorm))
    CASE DEFAULT
      STOP 'RhoFluxVar for GVEC eq. not between 0 or 3'
    END SELECT
    
    MHDEQdata(  1,iNode)=Density
    MHDEQdata(  2,iNode)=MHDEQdata(2,iNode)*mu0 !pressure transformed to mu0=1
!    MHDEQdata( 3:5,iNode)=Bcart(:)
!    MHDEQdata(   6,iNode)=chi_int !poloidal flux
!    MHDEQdata(   7,iNode)=Phi_int !toroidal flux
!    MHDEQdata(8:10,iNode)=Acart(:)

  END DO !iNode=1,nTotal
  CALL ProgressBar(maxNode,nTotal)
END DO !istep
!WRITE(UNIT_stdOut,'(A,4(1X,E12.5))')'  Rmin/Zmax, Zmin/Zmax ', Rmin,Rmax,Zmin,Zmax 
WRITE(UNIT_stdOut,'(A,2(1X,E12.5))')'  xmin/xmax', MINVAL(x_out(1,:)),MAXVAL(x_out(1,:))
WRITE(UNIT_stdOut,'(A,2(1X,E12.5))')'  ymin/ymax', MINVAL(x_out(2,:)),MAXVAL(x_out(2,:))
WRITE(UNIT_stdOut,'(A,2(1X,E12.5))')'  zmin/zmax', MINVAL(x_out(3,:)),MAXVAL(x_out(3,:))

WRITE(UNIT_stdOut,'(A)')'  ...DONE.                             '

END SUBROUTINE MapToGVEC 


!===================================================================================================================================
!> Finalize GVEC module
!!
!===================================================================================================================================
SUBROUTINE FinalizeGVEC 
! MODULES
USE MOD_GVEC_Vars
#ifdef PP_GVEC
USE MODgvec_gvec_to_hopr,ONLY:Finalize_gvec_to_hopr
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#ifdef PP_GVEC
  CALL Finalize_gvec_to_hopr()
#endif

END SUBROUTINE FinalizeGVEC

END MODULE MOD_GVEC

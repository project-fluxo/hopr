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
MODULE MOD_MHDEQ
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

INTERFACE InitMHDEQ 
  MODULE PROCEDURE InitMHDEQ
END INTERFACE

! allow different dimensions of input/output arrays
!INTERFACE MapToMHDEQ
!  MODULE PROCEDURE MapToMHDEQ
!END INTERFACE

INTERFACE FinalizeMHDEQ 
  MODULE PROCEDURE FinalizeMHDEQ
END INTERFACE

PUBLIC::InitMHDEQ
PUBLIC::MapToMHDEQ
PUBLIC::FinalizeMHDEQ
!===================================================================================================================================

CONTAINS
SUBROUTINE InitMHDEQ 
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,abort
USE MOD_ReadInTools,ONLY:GETINT,GETREALARRAY
USE MOD_MHDEQ_Vars
USE MOD_VMEC, ONLY:InitVMEC
USE MOD_GVEC, ONLY:InitGVEC
USE MOD_Solov, ONLY:InitSolov
USE MOD_cyl1d, ONLY:Initcyl1d
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
WRITE(UNIT_stdOut,'(A)')'INIT MHD EQUILIBRIUM INPUT ...'
whichEquilibrium    = GETINT('whichEquilibrium','0')   
IF(WhichEquilibrium.EQ.0) THEN 
  WRITE(UNIT_stdOut,'(A)')'... NOTHING TO BE DONE'
  RETURN
END IF
SELECT CASE(whichEquilibrium)
CASE(1)
  useMHDEQ=.TRUE.
  nVarMHDEQ =10
  WRITE(*,*)'Using VMEC as equilibrium solution...'
  CALL InitVMEC()
CASE(2)
  useMHDEQ=.TRUE.
  nVarMHDEQ =10
  WRITE(*,*)'Using Soloviev as equilibrium solution...'
  CALL InitSolov()
CASE(3)
  useMHDEQ=.TRUE.
  nVarMHDEQ =10
  WRITE(*,*)'Using GVEC as equilibrium solution...'
  CALL InitGVEC()
CASE(4)
  useMHDEQ=.TRUE.
  nVarMHDEQ =13
  WRITE(*,*)'Using CYL1D as equilibrium solution...'
  CALL InitCyl1d()
CASE DEFAULT
  WRITE(*,*)'WARNING: No Equilibrium solution for which Equilibrium= ', whichEquilibrium
  STOP
END SELECT
  !density coefficients of the polynomial coefficients: rho_1+rho_2*x + rho_3*x^2 ...
  nRhoCoefs=GETINT("nRhoCoefs","1")
  ALLOCATE(RhoCoefs(nRhoCoefs))
  IF(nRhoCoefs.EQ.1)THEN !default = 1.
    RhoFluxVar=GETINT("RhoFluxVar","0") ! dependant variable: =0: psinorm (tor. flux), =1:chinorm (pol. flux)
    RhoCoefs=GETREALARRAY("RhoCoefs",nRhoCoefs,"1.0")
  ELSE
    RhoFluxVar=GETINT("RhoFluxVar") ! dependant variable: =0: psinorm (tor. flux), =1:chinorm (pol. flux)
    RhoCoefs=GETREALARRAY("RhoCoefs",nRhoCoefs)
  END IF
  InputCoordSys=GETINT("MHDEQ_inputCoordSys","0")
  ! =0: x_in(1:3) are (x,y,z) coordinates in a cylinder of size r=[0;1], z=[0;1]
  ! =1: x_in(1:3) are (r,zeta,theta) coordinates r= [0;1], zeta= [0;1], theta=[0;1]

WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE InitMHDEQ


SUBROUTINE MapToMHDEQ(nTotal,x_in,x_out,MHDEQdata)
!===================================================================================================================================
! Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from VMEC data. 
! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars, ONLY: nVarMHDEQ,MHDEQvarNames,whichEquilibrium,InputCoordSys
USE MOD_VMEC, ONLY:MapToVMEC
USE MOD_GVEC, ONLY:MapToGVEC
USE MOD_Solov, ONLY:MapToSolov
USE MOD_Cyl1d, ONLY:MapToCyl1d
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nTotal         ! total number of points
REAL, INTENT(IN)   :: x_in(3,nTotal) ! input coordinates represent a cylinder: 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: x_out(3,nTotal) ! mapped x,y,z coordinates with vmec data
REAL,INTENT(OUT)   :: MHDEQdata(nVarMHDEQ,nTotal) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!===================================================================================================================================
SELECT CASE(whichEquilibrium)
CASE(1)
  CALL MapToVMEC(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
CASE(2)
  CALL MapToSolov(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
CASE(3)
  CALL MapToGVEC(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
CASE(4)
  CALL MapToCyl1d(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
END SELECT
WRITE(*,'(A21,2(1X,A21))')'MHDeq Var ',' Min' ,' Max'
DO i=1,nVarMHDEQ
  WRITE(*,'(A21,2(1X,E21.11))')MHDEQvarNames(i),MINVAL(MHDEQdata(i,:)),MAXVAL(MHDEQdata(i,:))
END DO
END SUBROUTINE MapToMHDEQ 

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeMHDEQ 
! MODULES
USE MOD_MHDEQ_Vars
USE MOD_VMEC,  ONLY:FinalizeVMEC
USE MOD_GVEC,  ONLY:FinalizeGVEC
USE MOD_Solov, ONLY:FinalizeSolov

IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(whichEquilibrium)
CASE(1)
  CALL FinalizeVMEC()
CASE(2)
  CALL FinalizeSolov()
CASE(3)
  CALL FinalizeGVEC()
END SELECT

DEALLOCATE(MHDEQoutdataGL)
DEALLOCATE(MHDEQdataEq)
DEALLOCATE(RhoCoefs)


END SUBROUTINE FinalizeMHDEQ

END MODULE MOD_MHDEQ

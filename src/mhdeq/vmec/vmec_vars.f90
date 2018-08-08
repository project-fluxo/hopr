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
MODULE MOD_VMEC_Vars
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------

INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(15, 307) !working precision 

! GLOBAL VARIABLES 
LOGICAL                 :: useVMEC                   !< main switch
LOGICAL                 :: useSFL                    !< use straight-field line coordinates
CHARACTER(LEN = 256)    :: VMECdataFile
INTEGER,ALLOCATABLE     :: xmAbs(:)                  !< |xm(iMode)|, 1 for m=0, 2 for even, 3 for odd
REAL(wp),ALLOCATABLE    :: Phi_prof(:)               !< TOROIDAL flux profile (called phi in VMEC)
REAL(wp),ALLOCATABLE    :: Phinorm_prof(:)           !< normalized TOROIDAL flux profile 
REAL(wp),ALLOCATABLE    :: chi_prof(:)               !< POLOIDAL flux profile (called chi in VMEC)

REAL(wp),ALLOCATABLE    :: rho(:)                    !< := sqrt(phinorm) at all flux surface 
REAL(wp),ALLOCATABLE    :: pres_Spl(:,:)             !< Spline coefficients in (rho) for Pressure, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: Phi_Spl(:,:)              !< Spline coefficients in (rho) for Phi, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: chi_Spl(:,:)              !< Spline coefficients in (rho) for chi, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: iota_Spl(:,:)             !< Spline coefficients in (rho) for chi, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: Rmnc_Spl(:,:,:)           !< modified spline coefficients R cosine, (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: Rmns_Spl(:,:,:)           !< modified spline coefficients R sine,   (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: lmnc_Spl(:,:,:)           !< modified spline coefficients of lambda cosine , (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: lmns_Spl(:,:,:)           !< modified spline coefficients of lambda sine,   (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: Zmnc_Spl(:,:,:)           !< modified spline coefficients of Z cosine, (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: Zmns_Spl(:,:,:)           !< modified spline coefficients of Z sine,   (1:4,iFlux,iMode)

!===================================================================================================================================
END MODULE MOD_VMEC_Vars


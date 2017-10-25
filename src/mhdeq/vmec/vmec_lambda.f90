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
MODULE MOD_VMEC_Lambda

  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(15, 307) !working precision 
! ----------------------------------------------------------------------
!> recompute lambda at the full flux surfaces via an elliptic solution on 
!> each flux surface
!> integration over angles done wih trapezoidal rule, cosine and sine are
!> sampled at theta/zeta = [0,...,k/(2pi),...,(np-1)/(2pi)]
! ----------------------------------------------------------------------
INTERFACE RecomputeLambda
  MODULE PROCEDURE RecomputeLambda
END INTERFACE

CONTAINS

SUBROUTINE RecomputeLambda(np_m,np_n)
!=====================================================================
! INPUT/OUTPUT VARIABLES
USE MOD_Basis1D,    ONLY: SOLVE
USE MOD_VMEC_vars,  ONLY: xmabs,rho,Rmnc_Spl,Rmns_Spl,Zmnc_Spl,Zmns_Spl
USE MOD_VMEC_Readin,ONLY: Rmnc,Rmns,Zmnc,Zmns,lmnc,lmns,iotaf,lasym,xm,xn
USE MOD_VMEC_Readin,ONLY: nfp,nFluxVMEC,mn_mode,twoPi
USE SPLINE1_MOD, ONLY: SPLINE1_EVAL
!---------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER          :: np_m !< number of points for integration in theta 
  INTEGER          :: np_n !< number of points for integration in zeta
!---------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                                   :: iFlux,iMode,jMode,mn0Mode
  INTEGER                                   :: i_m,i_n,iGuess
  REAL(wp),DIMENSION(1:np_m,1:np_n,mn_mode) :: cosMN,sinMN
  REAL(wp),DIMENSION(1:np_m,1:np_n)         :: R,dRdrho,dRdtheta,dRdzeta
  REAL(wp),DIMENSION(1:np_m,1:np_n)         :: dZdrho,dZdtheta,dZdzeta
  REAL(wp),DIMENSION(1:np_m,1:np_n)         :: J,g_tt,g_tz,g_zz 
  REAL(wp),ALLOCATABLE                      :: Amat(:,:),RHS(:),lambda(:),Adiag(:)
  REAL(wp)                                  :: rho_p,rhom,drhom,splOut(3) !for weighted spline interpolation
!=====================================================================
  WRITE(*,'(4X,A)')'RECOMPUTE LAMBDA ON FULL GRID...'
  IF(lasym)THEN !cosine & sine
    ALLOCATE(Amat(2*mn_mode,2*mn_mode),RHS(2*mn_mode),lambda(2*mn_mode))
    !estimate of Adiag
  ELSE !only sine
    ALLOCATE(Amat(mn_mode,mn_mode),RHS(mn_mode),lambda(mn_mode))
  END IF
  !estimate of Adiag
  ALLOCATE(Adiag(1:mn_mode)) !Adiag( mn_mode+iMode)=Adiag(iMode)
  DO iMode=1,mn_mode
    Adiag(iMode)=MAX(1.0_wp,((REAL(nfp,wp)*xm(iMode))**2+(xn(iMode))**2) )*np_m*np_n
  END DO

  DO iMode=1,mn_mode
    IF((xm(iMode).EQ.0).AND.(xn(iMode).EQ.0)) mn0Mode=iMode
  END DO

  DO i_m=1,np_m; DO i_n=1,np_n
    cosMN(i_m,i_n,:)=COS(twoPi*( xm(:)*(REAL(i_m-1,wp)/REAL(np_m,wp)) & 
                                -xn(:)*(REAL(i_n-1,wp)/REAL(nfp*np_n,wp))))
    sinMN(i_m,i_n,:)=SIN(twoPi*( xm(:)*(REAL(i_m-1,wp)/REAL(np_m,wp)) &
                                -xn(:)*(REAL(i_n-1,wp)/REAL(nfp*np_n,wp))))
  END DO; END DO !i_m,i_n
  cosMN(:,:,mn0Mode)=1.0_wp
  sinMN(:,:,mn0Mode)=0.0_wp

  DO iFlux=2,nFluxVMEC
    WRITE(0,'(I6,A4,I6,A27,A1)',ADVANCE='NO')iFlux, ' of ',nFluxVMEC,' flux surfaces evaluated...',ACHAR(13)
    rho_p=rho(iFlux) !rho=sqrt(psi_norm)
    dRdrho   = 0.0_wp
    dRdrho   = 0.0_wp
    R        = 0.0_wp
    dRdtheta = 0.0_wp
    dRdzeta  = 0.0_wp
    dZdrho   = 0.0_wp
    dZdtheta = 0.0_wp
    dZdzeta  = 0.0_wp
    DO iMode=1,mn_mode
      SELECT CASE(xmabs(iMode))
      CASE(0)
        rhom=1.0_wp
        drhom=0.0_wp
      CASE(1)
        rhom=rho_p
        drhom=1.0_wp
      CASE(2)
        rhom=rho_p*rho_p
        drhom=2.0_wp*rho_p
      CASE DEFAULT
        rhom=rho_p**xmabs(iMode)
        drhom=REAL(xmabs(iMode) , wp)*rho_p**(xmabs(iMode)-1)
      END SELECT
      R(:,:)       =       R(:,:)          +Rmnc(iMode,iFlux)*cosMN(:,:,iMode)
      dRdtheta(:,:)=dRdtheta(:,:)-xm(iMode)*Rmnc(iMode,iFlux)*sinMN(:,:,iMode)
      dRdzeta(:,:) = dRdzeta(:,:)+xn(iMode)*Rmnc(iMode,iFlux)*sinMN(:,:,iMode)

      dZdtheta(:,:)=dZdtheta(:,:)+xm(iMode)*Zmns(iMode,iFlux)*cosMN(:,:,iMode)
      dZdzeta(:,:) = dZdzeta(:,:)-xn(iMode)*Zmns(iMode,iFlux)*cosMN(:,:,iMode)

      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Rmnc_Spl(:,:,iMode),iGuess,splout) 
      dRdrho(:,:)=dRdrho(:,:)+ (rhom*splout(2)+splout(1)*drhom)*CosMN(:,:,iMode)
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Zmns_Spl(:,:,iMode),iGuess,splout) 
      dZdrho(:,:)=dZdrho(:,:)+ (rhom*splout(2)+splout(1)*drhom)*SinMN(:,:,iMode)
      IF(lasym)THEN
        R(:,:)       =       R(:,:)+          Rmns(iMode,iFlux)*sinMN(:,:,iMode)
        dRdtheta(:,:)=dRdtheta(:,:)+xm(iMode)*Rmns(iMode,iFlux)*cosMN(:,:,iMode)
        dRdzeta(:,:) = dRdzeta(:,:)-xn(iMode)*Rmns(iMode,iFlux)*cosMN(:,:,iMode)

        dZdtheta(:,:)=dZdtheta(:,:)-xm(iMode)*Zmnc(iMode,iFlux)*sinMN(:,:,iMode)
        dZdzeta(:,:) = dZdzeta(:,:)+xn(iMode)*Zmnc(iMode,iFlux)*sinMN(:,:,iMode)

        CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Rmns_Spl(:,:,iMode),iGuess,splout) 
        dRdrho(:,:)=dRdrho(:,:)+ (rhom*splout(2)+splout(1)*drhom)*SinMN(:,:,iMode)

        CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Zmnc_Spl(:,:,iMode),iGuess,splout) 
        dZdrho(:,:)= dZdrho(:,:)+ (rhom*splout(2)+splout(1)*drhom)*CosMN(:,:,iMode)
      END IF !lasym
    END DO!iMode
    !J=sqrtG=R*(dRdtheta*dZds-dRds*dZdtheta) =  1/(2*rho_p)*R*(dRdtheta*dZdrho-dRdrho*dZdtheta) 
    J=R*(dRdtheta*dZdrho-dRdrho*dZdtheta) !1/(2rho_p) is constant on flux, not needed
    IF(MINVAL(ABS(J)).LT.1.0e-12) STOP 'Jacobian near zero! J< 1.0e-12'
    !account for 1/J here
    g_tt =(dRdtheta*dRdtheta+dZdtheta*dZdtheta      )/J
    g_tz =(dRdtheta*dRdzeta +dZdtheta*dZdzeta       )/J !=g_zt
    g_zz =(dRdzeta *dRdzeta +dZdzeta *dZdzeta + R**2)/J
    
    
    
    !lambda(1:mn_mode)=lmns(:,iFlux)
    DO jMode=1,mn_mode
      DO iMode=1,mn_mode
        ! 1/J ( (g_thet,zeta  dlambdaSIN_dzeta -g_zeta,zeta  dlambdaSIN_dthet) dsigmaSIN_dthet
        !      +(g_thet,thet  dlambdaSIN_dzeta -g_zeta,thet  dlambdaSIN_dthet) dsigmaSIN_dzeta)
        Amat(iMode,jMode) = &
!                          SUM(( ( g_tz(:,:)*(-xn(iMode)) -g_zz(:,:)*( xm(iMode)))*( xm(jMode))    &  
!                               +( g_tt(:,:)*(-xn(iMode)) -g_tz(:,:)*( xm(iMode)))*(-xn(jMode)) )  &  
!optimized
                          SUM((  g_tz(:,:)*(((-xn(iMode))*( xm(jMode))) -           (( xm(iMode))*(-xn(jMode))))   &  
                               + g_tt(:,:)* ((-xn(iMode))*(-xn(jMode))) - g_zz(:,:)*(( xm(iMode))*( xm(jMode))) )  &  
                                       *cosMN(:,:,iMode)*cosMN(:,:,jMode) )/Adiag(iMode)
      END DO !iMode
      ! 1/J( (iota g_thet,zeta  + g_zeta,zeta )  dsigmaSIN_dthet
      !     +(iota g_thet,thet  + g_zeta,thet )  dsigmaSIN_dzeta
      rhs(jMode)      =   &
!                         SUM(( ( g_tz(:,:)*iotaf(iFlux) +g_zz(:,:)  )*( xm(jMode))                &  
!                              +( g_tt(:,:)*iotaf(iFlux) +g_tz(:,:)  )*(-xn(jMode)))               &  
!optimized
                         SUM( ( g_tz(:,:)*(iotaf(iFlux)*( xm(jMode))  +           (-xn(jMode)))   &  
                               +g_tt(:,:)*(iotaf(iFlux)*(-xn(jMode))) + g_zz(:,:)*( xm(jMode)) )  & 
                                                  *cosMN(:,:,jMode) )/Adiag(jMode)
    END DO!jMode
    Amat(mn0mode,:)      =0.0_wp
    Amat(mn0mode,mn0mode)=1.0_wp
    RHS(         mn0mode)=0.0_wp
    IF(lasym)THEN
      !lambda(1:mn_mode)=lmns(:,iFlux)
      !lambda(mn_mode+1:2*mn_mode)=lmnc(:,iFlux)
      do jmode=1,mn_mode
        do imode=1,mn_mode
          ! 1/J ( (g_thet,thet dlambdaCOS_dzeta-g_thet,zeta  dlambdaCOS_dthet) dsigmaSIN_dthet
          !      +(g_zeta,thet dlambdaCOS_dzeta-g_zeta,zeta  dlambdaCOS_dthet) dsigmaSIN_dzeta)
          Amat(mn_mode+iMode,        jMode)= &
!                          SUM(( ( g_tz(:,:)*( xn(iMode)) -g_zz(:,:)*(-xm(iMode)))*( xm(jMode))    &  
!                               +( g_tt(:,:)*( xn(iMode)) -g_tz(:,:)*(-xm(iMode)))*(-xn(jMode)) )  &  
!optimized
                          SUM((  g_tz(:,:)*(( xn(iMode))*( xm(jMode))  -           (-xm(iMode))*(-xn(jMode)))    &  
                               + g_tt(:,:)*(( xn(iMode))*(-xn(jMode))) -g_zz(:,:)*((-xm(iMode))*( xm(jMode))) )  &  
                                      *sinMN(:,:,iMode)*cosMN(:,:,jMode) )/Adiag(iMode)
          ! 1/J ( (g_thet,thet dlambdaSIN_dzeta-g_thet,zeta  dlambdaSIN_dthet) dsigmaCOS_dthet
          !      +(g_zeta,thet dlambdaSIN_dzeta-g_zeta,zeta  dlambdaSIN_dthet) dsigmaCOS_dzeta)
          Amat(        iMode,mn_mode+jMode)= &
!                          SUM(( ( g_tz(:,:)*(-xn(iMode)) -g_zz(:,:)*( xm(iMode)))*(-xm(jMode))    &  
!                               +( g_tt(:,:)*(-xn(iMode)) -g_tz(:,:)*( xm(iMode)))*( xn(jMode)) )  &  
!optimized
                          SUM((  g_tz(:,:)*((-xn(iMode))*(-xm(jMode))  -           ( xm(iMode))*( xn(jMode)))    &  
                               + g_tt(:,:)*((-xn(iMode))*( xn(jMode))) -g_zz(:,:)*(( xm(iMode))*(-xm(jMode))))  &  
                                      *cosMN(:,:,iMode)*sinMN(:,:,jMode) )/Adiag(iMode)
          ! 1/J ( (g_thet,thet dlambdaCOS_dzeta-g_thet,zeta  dlambdaCOS_dthet) dsigmaCOS_dthet
          !      +(g_zeta,thet dlambdaCOS_dzeta-g_zeta,zeta  dlambdaCOS_dthet) dsigmaCOS_dzeta)
          Amat(mn_mode+iMode,mn_mode+jMode)= &
!                          SUM(( ( g_tz(:,:)*( xn(iMode)) -g_zz(:,:)*(-xm(iMode)))*(-xm(jMode))    &  
!                               +( g_tt(:,:)*( xn(iMode)) -g_tz(:,:)*(-xm(iMode)))*( xn(jMode)) )  &  
!optimized
                          SUM((  g_tz(:,:)*(( xn(iMode))*(-xm(jMode))  -           (-xm(iMode))*( xn(jMode)))    &  
                               + g_tt(:,:)*(( xn(iMode))*( xn(jMode))) -g_zz(:,:)*((-xm(iMode))*(-xm(jMode))) )  &  
                                      *sinMN(:,:,iMode)*sinMN(:,:,jMode) )/Adiag(iMode)


        END DO !iMode
        ! 1/J( (iota g_thet,thet+ g_thet,zeta)  dsigmaCOS_dthet
        !     +(iota g_zeta,thet+ g_zeta,zeta)  dsigmaCOS_dzeta
!        RHS(mn_mode+jMode) = SUM(( ( g_tz(:,:)*iotaf(iFlux) +g_zz(:,:)  )*(-xm(jMode))               &  
!                                  +( g_tt(:,:)*iotaf(iFlux) +g_tz(:,:)  )*( xn(jMode)))              &  
!optimized
        RHS(mn_mode+jMode) = SUM((   g_tz(:,:)*((iotaf(iFlux)*(-xm(jMode))) +          ( xn(jMode)))     &  
                                  +  g_tt(:,:)* (iotaf(iFlux)*( xn(jMode))) +g_zz(:,:)*(-xm(jMode))   )  &  
                                                                    *sinMN(:,:,jMode) )/Adiag(jMode)

      END DO!jMode
      Amat(mn_mode+mn0Mode,:              )=0.0_wp
      Amat(mn_mode+mn0Mode,mn_mode+mn0Mode)=1.0_wp
      RHS(                 mn_mode+mn0Mode)=0.0_wp
    END IF !lasym
   
    lambda=SOLVE(Amat,RHS)  
!    WRITE(*,'(A,I6,2E21.11)')'DEBUG',iFlux,MAXVAL(ABS(lambda(1:mn_mode)-lmns(:,iFlux)))
!    DO iMode=1,mn_mode
!      WRITE(*,'(A,I6,I6,3E21.11)')'DEBUG',NINT(xm(iMode)),NINT(xn(iMode)),lambda(iMode),lmns(iMode,iFlux),lambda(iMode)/(1.0e-13+lmns(iMode,iFlux))
!    END DO
!    WRITE(*,*)'***************************'

    !overwrite
    lmns(:,iFlux)=lambda(1:mn_mode)
    IF(lasym)THEN
!    WRITE(*,'(A,I6,2E21.11)')'DEBUG',iFlux,MAXVAL(ABS(lambda(mn_mode+1:2*mn_mode)-lmnc(:,iFlux)))
!    DO iMode=1,mn_mode
!      WRITE(*,'(A,I6,I6,3E21.11)')'DEBUG',NINT(xm(iMode)),NINT(xn(iMode)),lambda(mn_mode+iMode),lmnc(iMode,iFlux) &
!                            ,lambda(mn_mode+iMode)/(1.0e-13+lmnc(iMode,iFlux))
!    END DO
!    WRITE(*,*)'###########################'

      lmnc(:,iFlux)=lambda(mn_mode+1:2*mn_mode)
    END IF
  END DO !iFlux
  !on axis
  lmns(:,1)=0.0_wp
  IF(lasym)THEN
    lmnc(:,1)=0.0_wp
  END IF

  DEALLOCATE(lambda,Amat,RHS,Adiag)
  WRITE(*,'(A)')'                                                       '
  WRITE(*,'(4X,A)')'...DONE.'
END SUBROUTINE RecomputeLambda

END MODULE MOD_VMEC_lambda 

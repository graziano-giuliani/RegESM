!=======================================================================
! Regional Earth System Model (RegESM)
! Copyright (c) 2013-2017 Ufuk Turuncoglu
! Licensed under the MIT License.
!=======================================================================
#define FILENAME "mod_esmf_cop_void.F90"

!-----------------------------------------------------------------------
! COP gridded component code
!-----------------------------------------------------------------------

module mod_esmf_cop

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only : NUOPC_ModelGet,                             &
      NUOPC_SetServices            => SetServices,                    &
      NUOPC_Label_Advertise        => label_Advertise,                &
      NUOPC_Label_RealizeProvided  => label_RealizeProvided,          &
      NUOPC_Label_SetClock         => label_SetClock,                 &
      NUOPC_Label_Advance          => label_Advance,                  &
      NUOPC_Label_DataInitialize   => label_DataInitialize

  implicit none
  private

!-----------------------------------------------------------------------
! Public subroutines
!-----------------------------------------------------------------------

  public :: COP_SetServices

  contains

    subroutine COP_SetServices(gcomp, rc)
      implicit none
      type(ESMF_GridComp) :: gcomp
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      call NUOPC_CompDerive(gcomp, NUOPC_SetServices, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

      call NUOPC_CompAttributeSet(gcomp, name="Verbosity", value="high", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

      ! specialize model
      call NUOPC_CompSpecialize(gcomp, specLabel=NUOPC_label_Advertise, &
                                specRoutine=Advertise, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call NUOPC_CompSpecialize(gcomp, specLabel=NUOPC_label_RealizeProvided, &
                                specRoutine=Realize, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call NUOPC_CompSpecialize(gcomp, specLabel=NUOPC_label_SetClock, &
                                specRoutine=SetClock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call NUOPC_CompSpecialize(gcomp, specLabel=NUOPC_label_Advance, &
                                specRoutine=Advance, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

    end subroutine COP_SetServices

    subroutine Advertise(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc
      rc = ESMF_SUCCESS
    end subroutine Advertise

    subroutine Realize(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc
      rc = ESMF_SUCCESS
    end subroutine Realize

    subroutine SetClock(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc
      ! local variables
      type(ESMF_Clock)              :: clock
      type(ESMF_TimeInterval)       :: timeStep

      rc = ESMF_SUCCESS
      call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call ESMF_ClockSet(clock, timeStep=timeStep/8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
    end subroutine SetClock

    subroutine Advance(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc
      rc = ESMF_SUCCESS
    end subroutine Advance

end module mod_esmf_cop

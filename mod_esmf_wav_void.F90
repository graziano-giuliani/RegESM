!=======================================================================
! Regional Earth System Model (RegESM)
! Copyright (c) 2013-2019 Ufuk Turuncoglu
! Licensed under the MIT License.
!=======================================================================
#define FILENAME "mod_esmf_wav_void.F90"

!-----------------------------------------------------------------------
! WAV gridded component code
!-----------------------------------------------------------------------

module mod_esmf_wav
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

  public :: WAV_SetServices

  contains

    subroutine WAV_SetServices(model, rc)
      implicit none
      type(ESMF_GridComp) :: model
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      call NUOPC_CompDerive(model, NUOPC_SetServices, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

      call NUOPC_CompAttributeSet(model, name="Verbosity", value="high", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

      ! specialize model
      call NUOPC_CompSpecialize(model, specLabel=NUOPC_label_Advertise, &
                                specRoutine=Advertise, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call NUOPC_CompSpecialize(model, specLabel=NUOPC_label_RealizeProvided, &
                                specRoutine=Realize, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call NUOPC_CompSpecialize(model, specLabel=NUOPC_label_SetClock, &
                                specRoutine=SetClock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call NUOPC_CompSpecialize(model, specLabel=NUOPC_label_Advance, &
                                specRoutine=Advance, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

    end subroutine WAV_SetServices

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
      ! Test with half a day
      call ESMF_ClockSet(clock, timeStep=timeStep/2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
    end subroutine SetClock

    subroutine Advance(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc

      type(ESMF_Clock)     :: clock
      type(ESMF_State)     :: importState, exportState
      character(len=160)   :: msgString

      rc = ESMF_SUCCESS

      call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call ESMF_ClockPrint(clock, options="currTime", &
          preString="------>Advancing WAV from: ", unit=msgString, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call ESMF_ClockPrint(clock, options="stopTime", &
          preString="---------------------> to: ", unit=msgString, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
    end subroutine Advance

end module mod_esmf_wav

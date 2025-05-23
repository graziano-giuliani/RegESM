!=======================================================================
! Regional Earth System Model (RegESM)
! Copyright (c) 2013-2019 Ufuk Turuncoglu
! Licensed under the MIT License.
!=======================================================================
#define FILENAME "util/mod_types.F90"

!-----------------------------------------------------------------------
! Module for user defined types
!-----------------------------------------------------------------------

module mod_types

  use ESMF
  use NUOPC

  implicit none

!-----------------------------------------------------------------------
!   Constants
!     cf3 - 1/rhow (rhow is reference density of seawater in kg/m3)
!-----------------------------------------------------------------------

    real(ESMF_KIND_R8), parameter :: cp = 3985.0d0
    real(ESMF_KIND_R8), parameter :: rho0 = 1025.0d0
    real(ESMF_KIND_R8), parameter :: cf1 = rho0*cp
    real(ESMF_KIND_R8), parameter :: cf2 = 1.0d0/cf1
    real(ESMF_KIND_R8), parameter :: cf3 = 1.0d0/rho0
    real(ESMF_KIND_R8), parameter :: day2s = 1.0d0/86400.0d0
    real(ESMF_KIND_R8), parameter :: mm2m = 1.0d0/1000.0d0
    real(ESMF_KIND_R8), parameter :: pi = 4.0d0*atan(1.0d0)
    real(ESMF_KIND_R8), parameter :: pi2 = 2.0d0*pi
    real(ESMF_KIND_R8), parameter :: phi = 0.5d0*pi
    real(ESMF_KIND_R8), parameter :: D2R = PI/180.0d0
    real(ESMF_KIND_R8), parameter :: R2D = 1.0d0/D2R
    real(ESMF_KIND_R8), parameter :: RADIUS = 6371.0d3

    real(ESMF_KIND_R8), parameter :: MISSING_R8 = 1.0d20
    real(ESMF_KIND_R4), parameter :: MISSING_R4 = 1.0e20
    real(ESMF_KIND_R8), parameter :: TOL_R8 = MISSING_R8/2.0d0
    real(ESMF_KIND_R4), parameter :: TOL_R4 = MISSING_R4/2.0

    real(ESMF_KIND_I4), parameter :: ZERO_I4 = 0
    real(ESMF_KIND_R8), parameter :: ZERO_R4 = 0.0
    real(ESMF_KIND_R8), parameter :: ZERO_R8 = 0.0d0
    real(ESMF_KIND_R8), parameter :: ONE_R8 = 1.0d0

    integer(ESMF_KIND_I4), parameter :: MAPPED_MASK = 99
    integer(ESMF_KIND_I4), parameter :: UNMAPPED_MASK = 98

    integer, parameter :: MAX_MAPPED_GRID = 10000
!Laura per il Nord Adriatico, default MAX_MAPPED_GRID=1000

!-----------------------------------------------------------------------
!   RTM river point data type
!-----------------------------------------------------------------------

    type RTM_River
      logical :: asIndex
      integer :: isActive
      real(ESMF_KIND_R8)  :: eRadius
      real(ESMF_KIND_R8)  :: lat, lon
      integer :: iindex, jindex
      integer :: dir
      integer :: npoints
      real(ESMF_KIND_R8)  :: monfac(12)
      integer :: rootPet
      integer :: mapSize
      real(ESMF_KIND_R8)  :: mapArea
      real(ESMF_KIND_R8)  :: mapTable(3,MAX_MAPPED_GRID) ! i, j and weight
    end type RTM_River

!-----------------------------------------------------------------------
!   ESM generic field data type
!-----------------------------------------------------------------------

    type ESM_Field
      integer :: fid
      integer :: rank
      integer :: gtype
      integer :: itype
      character(len=100) :: short_name
      character(len=100) :: long_name
      character(len=100) :: units
      character(len=100) :: export_units
      real(ESMF_KIND_R8) :: scale_factor
      real(ESMF_KIND_R8) :: add_offset
      logical :: enable_integral_adj
      type(ESMF_RouteHandle) :: rhandle
    end type ESM_Field

!-----------------------------------------------------------------------
!   ESM generic mesh data type
!-----------------------------------------------------------------------

    type ESM_Mesh
      integer :: gid
      integer :: gtype
      real(ESMF_KIND_R8), allocatable :: glon(:,:)
      real(ESMF_KIND_R8), allocatable :: glat(:,:)
      real(ESMF_KIND_R8), allocatable :: gare(:,:)
      integer(ESMF_KIND_I4), allocatable :: gmsk(:,:)
    end type ESM_Mesh

!-----------------------------------------------------------------------
!   ESM high-level generic model data type
!-----------------------------------------------------------------------

    type ESM_Model
      character(len=100) :: name
      integer :: nPets
      logical :: modActive
      integer :: tile(2)
      integer :: haloWidth
      integer :: nLevs
      real(ESMF_KIND_R8), allocatable :: levs(:)
      integer, allocatable :: petList(:)
      integer :: isLand
      integer :: isOcean
      type(ESM_Mesh), allocatable :: mesh(:)
      type(ESM_Field), allocatable :: importField(:)
      type(ESM_Field), allocatable :: exportField(:)
      type(ESMF_Grid) :: grid
      type(ESMF_Grid) :: grid3d
    end type ESM_Model

!-----------------------------------------------------------------------
!   ESM high-level generic connector data type
!-----------------------------------------------------------------------

    type ESM_Conn
      character(len=100) :: name
      integer :: divDT
      integer :: nPets
      logical :: modActive
      logical :: modExtrapolation
      integer :: modInteraction
      integer, allocatable :: petList(:)
    end type ESM_Conn

!-----------------------------------------------------------------------
!   ESM grided component (models) holder
!-----------------------------------------------------------------------

    type(ESM_Model), allocatable, target :: models(:)

!-----------------------------------------------------------------------
!   ESM connector (coupler) holder
!-----------------------------------------------------------------------

    type(ESM_Conn), allocatable, target :: connectors(:,:)

!-----------------------------------------------------------------------
!   Number of gridded component or model
!-----------------------------------------------------------------------

    integer :: nModels

!-----------------------------------------------------------------------
!   ESM model indices
!-----------------------------------------------------------------------

    character(len=3) :: COMPDES(5) = (/'ATM','OCN','RTM','WAV','COP'/)
    integer, parameter :: Iatmos = 1
    integer, parameter :: Iocean = 2
    integer, parameter :: Iriver = 3
    integer, parameter :: Iwavee = 4
    integer, parameter :: Icopro = 5

!-----------------------------------------------------------------------
!   Interaction interfaces
!-----------------------------------------------------------------------

    character(len=3) :: IFACEDES(3) = (/'LND','OCN','ALL'/)
    integer, parameter :: Ioverlnd = 1
    integer, parameter :: Ioverocn = 2
    integer, parameter :: Ioverall = 3

!-----------------------------------------------------------------------
!     Staggered grid point indices
!     d --------- d   d --- v --- d
!     |           |   |           |
!     |     c     |   u     c     u
!     |           |   |           |
!     d --------- d   d --- v --- d
!     Arakawa - B     Arakawa - C
!     RegCM           ROMS (c = rho, d = psi)
!-----------------------------------------------------------------------

    character(len=6) :: GRIDDES(0:4) = &
        (/'N/A   ','CROSS ','DOT   ','U     ','V     '/)
    integer, parameter :: Inan    = 0
    integer, parameter :: Icross  = 1
    integer, parameter :: Idot    = 2
    integer, parameter :: Iupoint = 3
    integer, parameter :: Ivpoint = 4

!-----------------------------------------------------------------------
!   Interpolation type
!-----------------------------------------------------------------------

    character(len=4) :: INTPDES(0:4) = (/'NONE','BLIN','CONS','NS2D', 'ND2S'/)
    integer, parameter :: Inone  = 0
    integer, parameter :: Ibilin = 1
    integer, parameter :: Iconsv = 2
    integer, parameter :: Instod = 3
    integer, parameter :: Indtos = 4

!-----------------------------------------------------------------------
!   Running mode
!-----------------------------------------------------------------------

    character(len=10) :: RUNNDES(2) = (/'SEQUENTIAL','CONCURRENT'/)
    integer, parameter :: Iseq = 1
    integer, parameter :: Ipar = 2

!-----------------------------------------------------------------------
!   ESM connector (coupler) holder
!-----------------------------------------------------------------------

    type(RTM_River), allocatable, target :: rivers(:)

!-----------------------------------------------------------------------
!   ESM model parameters
!-----------------------------------------------------------------------

    character(ESMF_MAXSTR) :: config_fname="namelist.rc"
    character(ESMF_MAXSTR), allocatable :: coproc_fnames(:)
    character(ESMF_MAXSTR) :: petLayoutOption
    type(ESMF_Time) :: esmStartTime
    type(ESMF_Time) :: esmRestartTime
    type(ESMF_Time) :: esmStopTime
    type(ESMF_TimeInterval) :: esmTimeStep
    type(ESMF_Calendar) :: esmCal
    type(ESMF_Clock) :: esmClock

    integer :: runMod
    integer :: cplType
    integer :: debugLevel
    logical :: enablePerfCheck
    logical :: restarted
    integer :: riverOpt

    contains

!-----------------------------------------------------------------------
!   Find index of specified field
!-----------------------------------------------------------------------

    integer function get_varid(list, key)
      implicit none

      type(ESM_Field), intent(in) :: list(:)
      character(len=*) :: key

      integer :: i

      do i = 1, size(list, dim=1)
        if (trim(list(i)%short_name) == trim(key)) then
          get_varid = i
          return
        end if
      end do
      get_varid = -1
      return
    end function get_varid

end module mod_types

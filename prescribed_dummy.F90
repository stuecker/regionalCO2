!-------------------------------------------------------------------
! manages reading and interpolation of prescribed dummy variable
! search for dummy/DUMMY
! Based on the prescribed_ozone.F90 routine created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_dummy

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_dummy_init
  public :: prescribed_dummy_adv
  public :: write_prescribed_dummy_restart
  public :: read_prescribed_dummy_restart
  public :: has_prescribed_dummy
  public :: prescribed_dummy_register
  public :: init_prescribed_dummy_restart
  public :: prescribed_dummy_readnl

  logical :: has_prescribed_dummy = .false.
  character(len=8), parameter :: dummy_name = 'dummy'

  character(len=16)  :: fld_name = 'DUMMY'
  character(len=256) :: filename = ' '
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_dummy_register()
    use ppgrid,         only: pver, pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: oz_idx


! changed to read in 2D data field: 
      call pbuf_add_field(dummy_name,'physpkg',dtype_r8,(/pcols,1/),oz_idx)
! ORIGINAL:       call pbuf_add_field(dummy_name,'physpkg',dtype_r8,(/pcols,pver/),oz_idx)

  endsubroutine prescribed_dummy_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_dummy_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: ndx, istat
    character(len=32) :: specifier(1)
    
    if ( has_prescribed_dummy ) then
       if ( masterproc ) then
          write(iulog,*) 'dummy variable is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    specifier(1) = trim(dummy_name)//':'//trim(fld_name)


    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)


  end subroutine prescribed_dummy_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_dummy_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_dummy_readnl'

   character(len=16)  :: prescribed_dummy_name
   character(len=256) :: prescribed_dummy_file
   character(len=256) :: prescribed_dummy_filelist
   character(len=256) :: prescribed_dummy_datapath
   character(len=32)  :: prescribed_dummy_type
   logical            :: prescribed_dummy_rmfile
   integer            :: prescribed_dummy_cycle_yr
   integer            :: prescribed_dummy_fixed_ymd
   integer            :: prescribed_dummy_fixed_tod

   namelist /prescribed_dummy_nl/ &
      prescribed_dummy_name,      &
      prescribed_dummy_file,      &
      prescribed_dummy_filelist,  &
      prescribed_dummy_datapath,  &
      prescribed_dummy_type,      &
      prescribed_dummy_rmfile,    &
      prescribed_dummy_cycle_yr,  &
      prescribed_dummy_fixed_ymd, &
      prescribed_dummy_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_dummy_name     = fld_name
   prescribed_dummy_file     = filename
   prescribed_dummy_filelist = filelist
   prescribed_dummy_datapath = datapath
   prescribed_dummy_type     = data_type
   prescribed_dummy_rmfile   = rmv_file
   prescribed_dummy_cycle_yr = cycle_yr
   prescribed_dummy_fixed_ymd= fixed_ymd
   prescribed_dummy_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_dummy_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_dummy_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_dummy_name,     len(prescribed_dummy_name),     mpichar, 0, mpicom)
   call mpibcast(prescribed_dummy_file,     len(prescribed_dummy_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_dummy_filelist, len(prescribed_dummy_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_dummy_datapath, len(prescribed_dummy_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_dummy_type,     len(prescribed_dummy_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_dummy_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_dummy_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_dummy_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_dummy_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   fld_name   = prescribed_dummy_name
   filename   = prescribed_dummy_file
   filelist   = prescribed_dummy_filelist
   datapath   = prescribed_dummy_datapath
   data_type  = prescribed_dummy_type
   rmv_file   = prescribed_dummy_rmfile
   cycle_yr   = prescribed_dummy_cycle_yr
   fixed_ymd  = prescribed_dummy_fixed_ymd
   fixed_tod  = prescribed_dummy_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_dummy = .true.

end subroutine prescribed_dummy_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_dummy_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use cam_control_mod, only: aqua_planet
    use phys_control, only : cam_physpkg_is
! added mwco2 here
    use physconst,    only : mwdry, mwco2                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz                ! J/K/molecule
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk, pbuf_get_field, pbuf_set_field

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    integer :: c,ncol
    real(r8) :: to_mmr(pcols,pver)

    real(r8) :: molmass
    real(r8) :: amass
    real(r8) :: outdata(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)

    character(len=32) :: units_str

    if( .not. has_prescribed_dummy ) return

       molmass = mwco2
       amass   = mwdry

    call advance_trcdata( fields, file, state, pbuf2d )

    units_str = trim(to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))))

!$OMP PARALLEL DO PRIVATE (C, NCOL, OUTDATA, TO_MMR, TMPPTR, PBUF_CHNK)
    do c = begchunk,endchunk
       ncol = state(c)%ncol

       select case ( units_str )
       case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
          to_mmr(:ncol,:) = (molmass*1.e6_r8*boltz*state(c)%t(:ncol,:))/(amass*state(c)%pmiddry(:ncol,:))
       case ('kg/kg','mmr')
          to_mmr(:ncol,:) = 1._r8
       case ('mol/mol','mole/mole','vmr','fraction')
          to_mmr(:ncol,:) = molmass/amass
       case default
          write(iulog,*) 'prescribed_dummy_adv: units = ',trim(fields(1)%units) ,' are not recognized'
          call endrun('prescribed_dummy_adv: units are not recognized')
       end select

       pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
       call pbuf_get_field(pbuf_chnk, fields(1)%pbuf_ndx, tmpptr )

       tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)

       outdata(:ncol,:) = (amass/molmass)* tmpptr(:ncol,:) ! vmr MFS
    enddo

  end subroutine prescribed_dummy_adv

!-------------------------------------------------------------------

  subroutine init_prescribed_dummy_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_dummy', piofile, file )

  end subroutine init_prescribed_dummy_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_dummy_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_dummy_restart

!-------------------------------------------------------------------
  subroutine read_prescribed_dummy_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile
    
    call read_trc_restart( 'prescribed_dummy', piofile, file )

  end subroutine read_prescribed_dummy_restart

end module prescribed_dummy

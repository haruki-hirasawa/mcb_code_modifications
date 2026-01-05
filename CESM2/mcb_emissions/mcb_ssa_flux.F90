module mcb_ssa_flux

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use shr_const_mod,        only : pi => shr_const_pi
  ! for reading netCDF
  use pio,              only : pio_inq_dimid, pio_inq_dimlen, pio_get_var, &
                               file_desc_t, var_desc_t, pio_inq_vardimid, pio_inq_varndims, pio_nowrite, &
                               pio_inq_varid, pio_closefile

  use cam_pio_utils,    only : cam_pio_openfile
  use cam_logfile,      only : iulog
  use sslt_sections,    only : nsections
  use spmd_utils,           only : masterproc
  implicit none

  private

  save

  character(len=256) :: mcb_filename = 'mcb_filename'
  
  real(r8)              :: mcb_emis_scale
  real(r8)              :: days(12)
  real(r8), allocatable :: &
       mcb_ssa(:,:,:)

  !!! should add up to 1
  real(r8), allocatable :: mcb_dist(:)
       
  public :: read_flux_mcb, get_fi_mcb, mcb_readnl
  
  contains

!!!
! Read namelist for SSa emission file
!!!
  subroutine mcb_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'mcb_readnl'

      namelist /mcb_nl/ mcb_filename, mcb_emis_scale
      !-----------------------------------------------------------------------------

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'mcb_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, mcb_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

      ! Broadcast namelist variables
      call mpibcast(mcb_filename, len(mcb_filename), mpichar, 0, mpicom)
      call mpibcast(mcb_emis_scale, 1,                            mpir8,   0, mpicom)
      print *, mcb_emis_scale
   end subroutine mcb_readnl


!-----------------------------------------------------------------------
! Initialize MCB emissions
!-----------------------------------------------------------------------
  subroutine read_flux_mcb
  ! mmm
    use spmd_utils,    only : masterproc
    use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, &
         interp_type ! for interpolating to model grid
    use ioFileMod,     only : getfil
    !use mo_chem_utls,  only : get_spc_ndx, get_extfrc_ndx
    use phys_grid,     only : get_ncols_p, get_rlat_all_p, get_rlon_all_p, ngcols_p
    use ppgrid,        only : begchunk, endchunk, pcols ! info on physics grid
    use time_manager, only : get_calday
    use sslt_sections, only: nsections
    
    implicit none
    
    !-------------------------------------------------------------------
    ! local variables
    real(r8), parameter :: d2r=pi/180._r8, zero=0._r8, twopi=pi*2._r8
    integer  :: k, j, n !iterables
    integer  :: nlat, nlon, ntime, ndims, nssbin !dimensions
    integer  :: ierr ! error code
    
    type(file_desc_t) :: piofile ! file object
    type(var_desc_t) :: vid ! variable id
    integer  :: dimid_lat, dimid_lon, dimid_time, dimid_ssbin ! dimension ids
    integer  :: dimid(3)
    real(r8), allocatable :: lat(:), lon(:) ! dimensions from input array
    integer, parameter :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
         716, 816, 915, 1016, 1115, 1216 /)
    real(r8), allocatable :: mcb_ssa_in(:,:,:) ! input array from netCDF (assume 3d - time, lat, lon)
    real(r8) :: total(2), tmp ! for printing global mean emis
    real(r8) :: factor ! for printing global mean emis
    character(len=256) :: locfn !for netCDF read
    type(interp_type) :: lon_wgts, lat_wgts ! weights for interpolation
    real(r8) :: to_lats(pcols), to_lons(pcols) ! target lon/lat in physics grid space
    integer :: ncols, c ! columns in a given chunk   
    
   !-----------------------------------------------------------------------
    !	... Open NetCDF file
    !-----------------------------------------------------------------------
    call getfil (mcb_filename, locfn, 0)
    call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
    
    !-----------------------------------------------------------------------
    !       ... get time dimension
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'mon', dimid_time )
    ierr = pio_inq_dimlen( piofile, dimid_time, ntime )
    if( ntime /= 12 )then
       write(iulog,*) 'read_flux_mcb: number of months = ',ntime,'; expecting 12'
       call endrun
    end if
    
    !-----------------------------------------------------------------------
    !       ... get sea salt bins dimension
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'ssbin', dimid_ssbin )
    ierr = pio_inq_dimlen( piofile, dimid_ssbin, nssbin )
    if( ierr /= 0 ) then
       write(iulog,*) 'read_flux_mcb: ssbin allocation error = ',ierr
       call endrun
    end if
    
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
    ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'read_flux_mcb: lat allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'lat', vid )
    ierr = pio_get_var( piofile, vid, lat )
    lat(:nlat) = lat(:nlat) * d2r
    !-----------------------------------------------------------------------
    !       ... get longitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
    ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
    allocate( lon(nlon), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'read_flux_mcb: lon allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'lon', vid )
    ierr = pio_get_var( piofile, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r
    !-----------------------------------------------------------------------
    !       ... Set up regridding
    !-----------------------------------------------------------------------

    allocate( mcb_ssa_in(nlon,nlat,ntime), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'read_flux_mcb: mcb_ssa_in allocation error = ',ierr
       call endrun
    end if

    allocate( mcb_ssa(pcols,begchunk:endchunk,ntime), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'read_flux_mcb: mcb_ssa allocation error = ',ierr
       call endrun
    end if

    allocate( mcb_dist(nssbin), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'read_flux_mcb: mcb_ssa_dist allocation error = ',ierr
       call endrun
    end if
    
    !-----------------------------------------------------------------------
    !	... Read emissions
    !-----------------------------------------------------------------------
    ierr = pio_inq_varid( piofile, 'mcb_ssa', vid )
    ierr = pio_inq_varndims( piofile, vid, ndims )
    
    if( ndims /= 3 ) then
       write(iulog,*) 'read_flux_mcb: variable mcb_ssa has ndims = ',ndims,', expecting 3'
       call endrun
    end if
    
    ierr = pio_inq_vardimid( piofile, vid, dimid )
    
    if( dimid(1) /= dimid_lon  .or. dimid(2) /= dimid_lat .or.  dimid(3) /= dimid_time ) then
       write(iulog,*) 'read_flux_mcb: Dimensions in wrong order for variable mcb_ssa'
       write(iulog,*) '...      Expecting (lon, lat, time)'
       call endrun
    end if
    
    ! reading data
    ierr = pio_get_var( piofile, vid, &
         (/ 1, 1, 1/), &                    ! start
         (/ nlon, nlat, ntime /), &   ! count
         mcb_ssa_in )
         
    !-----------------------------------------------------------------------
    !	... Read size distribution
    !-----------------------------------------------------------------------
    ierr = pio_inq_varid( piofile, 'mcb_dist', vid )
    ierr = pio_inq_varndims( piofile, vid, ndims )
    
    if( ndims /= 1 ) then
       write(iulog,*) 'read_flux_mcb: variable mcb_dist has ndims = ',ndims,', expecting 1'
       call endrun
    end if
    
    ierr = pio_inq_vardimid( piofile, vid, dimid )
    
    if( dimid(1) /= dimid_ssbin ) then
       write(iulog,*) 'read_flux_mcb: Dimensions in wrong order for variable mcb_dist'
       write(iulog,*) '...      Expecting (ssbin)'
       call endrun
    end if
    
    ! reading data
    ierr = pio_get_var( piofile, vid, &
         (/ 1 /), &                    ! start
         (/ nssbin /), &   ! count
         mcb_dist )

     ! finished reading
     call pio_closefile( piofile )
     
    !-----------------------------------------------------------------------
    !	... Regrid emissions
    !-----------------------------------------------------------------------
    allocate( mcb_ssa(pcols,begchunk:endchunk,ntime), stat=ierr )

    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)

       do k = 1,ntime
          call lininterp(mcb_ssa_in(:,:,k), nlon, nlat, mcb_ssa(:,c,k), ncols, lon_wgts, lat_wgts)
       enddo
       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
    enddo

    deallocate( mcb_ssa_in, lon, lat, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'read_flux_mcb: Failed to deallocate mcb_ssa_in; ierr = ',ierr
       call endrun
    end if
    
    !--------------------------------------------------------
    ! ... initialize the monthly day of year times
    !--------------------------------------------------------

    do n = 1,12
       days(n) = get_calday( dates(n), 0 )
    end do
    if (masterproc) then
       write(iulog,*) 'mcb_init : days'
       write(iulog,'(1p,5g15.8)') days(:)
    endif

    !!! At this point, the file has been read
    
  end subroutine read_flux_mcb
  
  
  !! Returns the MCB flux
  subroutine get_fi_mcb (lchnk, ncol, flux)
    use time_manager, only : get_curr_calday

    implicit none
    integer, intent(in) :: lchnk, ncol
    real(r8), intent(inout) :: flux(:,:)
    
    ! local variables
    integer :: ssbin, i, m
    
    real(r8)      :: calday                   ! day of year including fraction
    real(r8)      :: dels
    integer       :: last
    integer       :: next

    ! Determine the calendar day.
    calday = get_curr_calday()

    !--------------------------------------------------------
    ! ... setup the time interpolation
    !--------------------------------------------------------
    if( calday < days(1) ) then
        next = 1
        last = 12
        dels = (365._r8 + calday - days(12)) / (365._r8 + days(1) - days(12))
        else if( calday >= days(12) ) then
        next = 1
        last = 12
        dels = (calday - days(12)) / (365._r8 + days(1) - days(12))
        else
        do m = 11,1,-1
           if( calday >= days(m) ) then
              exit
           end if
        end do
        last = m
        next = m + 1
        dels = (calday - days(m)) / (days(m+1) - days(m))
    end if

    dels = max( min( 1._r8,dels ),0._r8 )
    
    !! Need to do two things: sample at the correct month
    !! and scale by the emission distribution
    i = 1
    ssbin = 1
    do i = 1, ncol
        do ssbin = 1,nsections
            flux(i,ssbin) = ( mcb_ssa(i,lchnk,last) &
                            + dels * (mcb_ssa(i,lchnk,next) - mcb_ssa(i,lchnk,last)) ) &
                            * mcb_dist(ssbin) * mcb_emis_scale
            !if (flux(i,ssbin) > 0.0_r8) then
            !    print *, flux(i,ssbin)
            !end if
        end do
    end do
  end subroutine get_fi_mcb
  
end module mcb_ssa_flux


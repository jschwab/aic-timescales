! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

  use star_lib
  use star_def
  use const_def
  use crlibm_lib

  implicit none

  integer, parameter :: extra_info_alloc = 1
  integer, parameter :: extra_info_get = 2
  integer, parameter :: extra_info_put = 3

  ! these routines are called by the standard run_star check_model
contains

  subroutine how_many_other_mesh_fcns(id, n)
    integer, intent(in) :: id
    integer, intent(out) :: n
    type (star_info), pointer :: s
    integer :: ierr

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    n = 0
    if (s% x_logical_ctrl(1)) n = n+1
    if (s% x_logical_ctrl(2)) n = n+1
  end subroutine how_many_other_mesh_fcns

  function urca(s,k) result(res)
    type (star_info), pointer :: s
    integer, intent(in) :: k
    real(dp) :: res
    res = s% eta(k) / s% x_ctrl(2)
  end function urca

  subroutine urca_other_mesh_fcn_data( &
       id, nfcns, names, gval_is_xa_function, vals1, ierr)
    integer, intent(in) :: id
    integer, intent(in) :: nfcns
    character (len=*) :: names(:)
    logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
    real(dp), pointer :: vals1(:) ! =(nz, nfcns)
    integer, intent(out) :: ierr
    integer :: nz, k
    real(dp), pointer :: vals(:,:)
    real(dp) :: weight
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    nz = s% nz
    vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)

    ! ensure N zones per decade in radius
    if (s% x_logical_ctrl(1)) then
       names(1) = 'logR'
       gval_is_xa_function(1) = .false.
       do k=1,nz
          vals(k,1) = s% x_ctrl(1) * log10_cr(s% r(k))
       end do
    end if

    ! remesh based on chemical potential
    if (s% x_logical_ctrl(2)) then
       names(2) = 'urca1'
       gval_is_xa_function(1) = .false.
       do k=1,nz
          vals(k,2) = urca(s, k)
       end do
    end if

  end subroutine urca_other_mesh_fcn_data

  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
    s% how_many_other_mesh_fcns => how_many_other_mesh_fcns
    s% other_mesh_fcn_data => urca_other_mesh_fcn_data

    ! Uncomment these lines if you wish to use the functions in this file,
    ! otherwise we use a null_ version which does nothing.
    s% extras_startup => extras_startup
    s% extras_start_step => extras_start_step
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns  

    s% how_many_extra_history_header_items => how_many_extra_history_header_items
    s% data_for_extra_history_header_items => data_for_extra_history_header_items
    s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
    s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

    ! Once you have set the function pointers you want,
    ! then uncomment this (or set it in your star_job inlist)
    ! to disable the printed warning message,
    s% job% warn_run_star_extras =.false.

  end subroutine extras_controls


  integer function extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_startup = 0
    if (.not. restart) then
       call alloc_extra_info(s)
    else ! it is a restart
       call unpack_extra_info(s)
    end if
  end function extras_startup

  integer function extras_start_step(id, id_extra)
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_start_step = 0
  end function extras_start_step



  ! returns either keep_going, retry, backup, or terminate.
  integer function extras_check_model(id, id_extra)

    use chem_def

    integer, intent(in) :: id, id_extra
    integer :: i, k, i_burn_max, nz, ierr
    type (star_info), pointer :: s
    real(dp) :: max_diff_eta, orunaway
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going

    ! apply varcontrol-esque timestep limit to eta
    if (s% x_logical_ctrl(3)) then

       i = maxloc(abs(s% eta(1:s% nz) - s% eta_start(1:s% nz)), dim=1)
       max_diff_eta = s% eta(i) - s% eta_start(i)

       if (abs(max_diff_eta) .gt. s% x_ctrl(3)) then
          extras_check_model = redo
          s% dt = 0.5d0 * s% dt
          write(*,*) 'redo for delta eta limit'
       endif

    endif

    ! terminate at oxygen runaway
    if (s% x_logical_ctrl(4)) then
       k = s% max_eps_nuc_k
       orunaway = s% eps_nuc_categories(ioo,k) - s% non_nuc_neu(k)

       if ((orunaway .gt. 0) .and. (s%T(k) .gt. 8e8)) then
          extras_check_model = terminate
          s% termination_code = t_xtra2
          termination_code_str(t_xtra2) = 'O+O runaway'
          return
       end if

    end if

    ! if you want to check multiple conditions, it can be useful
    ! to set a different termination code depending on which
    ! condition was triggered.  MESA provides 9 customizeable
    ! termination codes, named t_xtra1 .. t_xtra9.  You can
    ! customize the messages that will be printed upon exit by
    ! setting the corresponding termination_code_str value.
    ! termination_code_str(t_xtra1) = 'my termination condition'

    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id, id_extra)
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns = 4
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)

    use chem_def

    integer, intent(in) :: id, id_extra, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k, nz

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    !note: do NOT add the extras names to history_columns.list
    ! the history_columns.list is only for the built-in log column options.
    ! it must not include the new column names you are adding here.

    nz = s% nz

    names = 'unset'
    vals = 0

    names(1) = 't_accretion'
    vals(1) = abs(s% star_mass / s% star_mdot) * secyer

    names(2) = 'center_t_heat'
    vals(2) =  s% cp(nz) * s% T(nz) / s% eps_nuc(nz)

    names(3) = 'center_t_compress'
    vals(3) = 1d0 / s% dlnd_dt(nz)

    names(4) = 'center_t_capture'
    vals(4) = s% dt / ((s% ye_start(nz) - s% ye(nz)) / s% ye(nz))

  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id, id_extra)
    use star_def, only: star_info
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 6
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
    use chem_def
    use eos_def
    use eos_lib
    use net_lib, only: net_work_size
    use star_def, only: star_info, maxlen_profile_column_name
    integer, intent(in) :: id, id_extra, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: i, j, k, ci, op_err, net_lwork
    logical :: okay

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    net_lwork = net_work_size(s% net_handle, ierr)

    !note: do NOT add the extra names to profile_columns.list
    ! the profile_columns.list is only for the built-in profile column options.
    ! it must not include the new column names you are adding here.

    ! here is an example for adding a profile column
    !if (n /= 1) stop 'data_for_extra_profile_columns'

    names(1) = 'nz_minus_k'
    do k = 1, nz
       vals(k,1) = (s% nz + 1) - k
    end do

    names(2) = 'eps_nuc_mc2'

    vals(:,2) = 0
    do k=1, nz
       do i=1, s% species
          ci = s% chem_id(i)
          ! dxdt(i) = chem_isos% Z_plus_N(ci)*dydt(i_rate, i)
          vals(k,2) = vals(k,2) - chem_isos% mass_excess(ci)*s% dxdt_nuc(i,k)/chem_isos% Z_plus_N(ci)
       end do
    end do

    vals(:,2) = vals(:,2) * Qconv

    ! use reported neu value
    names(3) = 'eps_nuc_neu'
    vals(1:nz,3) = s% eps_nuc_neu_total(1:nz)


    ! calculate eps_eos as a remainder
    names(4) = 'eps_nuc_eos'
    do k=1, nz
       vals(k,4) = s% eps_nuc(k) - (vals(k,2) - vals(k,3))
    end do

    names(5) = 'log_rate_r_n20_wk_f20'
    names(6) = 'log_rate_r1616'

    okay = .true.
!$OMP PARALLEL DO PRIVATE(k,op_err) 
    do k = 1, s% nz
       if (.not. okay) cycle
       op_err = 0
       call do1_net( &
            s, k, s% species, &
            s% num_reactions, net_lwork, &
            n, nz, vals, op_err)
       if (op_err /= 0) okay = .false.
    end do
!$OMP END PARALLEL DO        


  end subroutine data_for_extra_profile_columns


  subroutine do1_net( &
       s, k, species, num_reactions, net_lwork, &
       n, nz, vals, ierr)
    use rates_def, only: std_reaction_Qs, std_reaction_neuQs, i_rate
    use net_def, only: Net_Info
    use net_lib, only: net_get
    use chem_def, only: chem_isos, category_name
    use eos_def, only : i_eta
    use utils_lib,only: &
         is_bad_num, realloc_double, realloc_double3
    type (star_info), pointer :: s         
    integer, intent(in) :: k, species, num_reactions, net_lwork, n, nz
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr

    integer :: i, j, screening_mode
    real(dp) :: log10_rho, log10_T, alfa, beta, &
         d_eps_nuc_dRho, d_eps_nuc_dT, cat_factor
    real(dp), target :: net_work_ary(net_lwork)
    real(dp), pointer :: net_work(:)
    type (Net_Info), target :: net_info_target
    type (Net_Info), pointer :: netinfo

    character (len=100) :: message
    real(dp), pointer :: reaction_neuQs(:)
    integer :: sz
    real(dp) :: eps_nuc_factor

    logical, parameter :: dbg = .false.

    include 'formats'

    ierr = 0

    net_work => net_work_ary
    netinfo => net_info_target

    log10_rho = s% lnd(k)/ln10
    log10_T = s% lnT(k)/ln10

    screening_mode = get_screening_mode(s,ierr)         
    if (ierr /= 0) then
       write(*,*) 'unknown string for screening_mode: ' // trim(s% screening_mode)
       stop 'do1_net'
       return
    end if

    call net_get( &
         s% net_handle, .false., netinfo, species, num_reactions, s% xa(1:species,k), &
         s% T(k), log10_T, s% rho(k), log10_Rho, &
         s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), &
         s% eta(k), s% d_eos_dlnd(i_eta,k), s% d_eos_dlnT(i_eta,k), &
         s% rate_factors, s% weak_rate_factor, &
         std_reaction_Qs, std_reaction_neuQs, .false., .false., &
         s% eps_nuc(k), d_eps_nuc_dRho, d_eps_nuc_dT, s% d_epsnuc_dx(:,k), & 
         s% dxdt_nuc(:,k), s% dxdt_dRho(:,k), s% dxdt_dT(:,k), s% d_dxdt_dx(:,:,k), &
         screening_mode, s% theta_e(k), &
         s% eps_nuc_categories(:,k), &
         s% eps_nuc_neu_total(k), net_lwork, net_work, ierr)

    if (ierr /= 0) then
       write(*,*) 'do1_net: net_get failure for cell ', k
       return
    end if

    call show_stuff(s,k,net_lwork,net_work,n,nz,vals)

  end subroutine do1_net


  integer function get_screening_mode(s,ierr)
    use rates_lib, only: screening_option
    type (star_info), pointer :: s 
    integer, intent(out) :: ierr
    include 'formats'
    ierr = 0
    if (s% screening_mode_value >= 0) then
       get_screening_mode = s% screening_mode_value
       return
    end if
    get_screening_mode = screening_option(s% screening_mode, ierr)
    if (ierr /= 0) return
    s% screening_mode_value = get_screening_mode
    !write(*,2) 'get_screening_mode ' // &
    !   trim(s% screening_mode), get_screening_mode
  end function get_screening_mode



  subroutine show_stuff(s,k,lwork,work,n,nz,vals)
    use chem_def
    use rates_def
    use rates_lib, only: rates_reaction_id
    use net_lib, only: get_reaction_id_table_ptr, get_net_rate_ptrs
    use crlibm_lib, only: log10_cr
    type (star_info), pointer :: s         
    integer, intent(in) :: k, lwork, n, nz
    real(dp), pointer :: work(:)
    real(dp) :: vals(nz,n)

    integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
    integer :: i, j, ierr, species, num_reactions, rate_id
    real(dp), pointer, dimension(:) :: &
         rate_screened, rate_screened_dT, rate_screened_dRho, &
         rate_raw, rate_raw_dT, rate_raw_dRho

    include 'formats'

    ierr = 0
    num_reactions = s% num_reactions

    call get_net_rate_ptrs(s% net_handle, &
         rate_screened, rate_screened_dT, rate_screened_dRho, &
         rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
         ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in get_net_rate_ptrs'
       stop 1
    end if

    call get_reaction_id_table_ptr(s% net_handle, reaction_id, ierr) 
    if (ierr /= 0) return

    rate_id = rates_reaction_id('r_ne20_wk_f20')
    if (rate_id <= 0) then
       write(*,*) 'failed to find reaction rate id -- not valid name?'
       vals(k,5) = 0
       return
    end if
    vals(k,5) = get_rate(rate_id)

    rate_id = rates_reaction_id('r_1616')
    if (rate_id <= 0) then
       write(*,*) 'failed to find reaction rate id -- not valid name?'
       vals(k,6) = 0
       return
    end if
    vals(k,6) = get_rate(rate_id)


  contains

    real(dp) function get_rate(id)
      integer, intent(in) :: id
      integer :: j
      include 'formats'

      get_rate = -99
      do j=1,num_reactions
         if (reaction_id(j) /= rate_id) cycle
         !write(*,3) 'screened rate ' // trim(reaction_Name(reaction_id(j))), &
         !   j, k, rate_screened(j)
         if (rate_screened(j) < 1d-20) then
            get_rate = -99
         else
            get_rate = log10_cr(rate_screened(j)) ! or rate_raw(j)
         end if
         return
      end do

      write(*,*) 'failed to find reaction rate id -- not in current net?'
      get_rate = -99

    end function get_rate

  end subroutine show_stuff


  
  subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
    integer, intent(in) :: id, id_extra
    integer, intent(out) :: num_cols
    num_cols=0
  end subroutine how_many_extra_history_header_items

  subroutine data_for_extra_history_header_items( &
       id, id_extra, num_extra_header_items, &
       extra_header_item_names, extra_header_item_vals, ierr)
    integer, intent(in) :: id, id_extra, num_extra_header_items
    character (len=*), pointer :: extra_header_item_names(:)
    real(dp), pointer :: extra_header_item_vals(:)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

    !here is an example for adding an extra history header item
    !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
    !extra_header_item_names(1) = 'mixing_length_alpha'
    !extra_header_item_vals(1) = s% mixing_length_alpha
  end subroutine data_for_extra_history_header_items


  subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
    integer, intent(in) :: id, id_extra
    integer, intent(out) :: num_cols
    num_cols = 0
  end subroutine how_many_extra_profile_header_items

  subroutine data_for_extra_profile_header_items( &
       id, id_extra, num_extra_header_items, &
       extra_header_item_names, extra_header_item_vals, ierr)
    integer, intent(in) :: id, id_extra, num_extra_header_items
    character (len=*), pointer :: extra_header_item_names(:)
    real(dp), pointer :: extra_header_item_vals(:)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

    !here is an example for adding an extra profile header item
    !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
    !extra_header_item_names(1) = 'mixing_length_alpha'
    !extra_header_item_vals(1) = s% mixing_length_alpha
  end subroutine data_for_extra_profile_header_items


  ! returns either keep_going or terminate.
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id, id_extra)
    type (star_info), pointer :: s
    integer, intent(in) :: id, id_extra
    integer :: ierr
    integer :: f
    extras_finish_step = keep_going

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    call store_extra_info(s)

    ! MESA provides a number of variables that make it easy to get user input.
    ! these are part of the star_info structure and are named
    ! x_character_ctrl, x_integer_ctrl, x_logical_ctrl, and x_ctrl.
    ! by default there are num_x_ctrls, which defaults to 100, of each.
    ! they can be specified in the controls section of your inlist.

    f = s% x_integer_ctrl(1)

    ! MESA also provides a number variables that are useful for implementing
    ! algorithms which require a state. if you just use these variables
    ! restarts, retries, and backups will work without doing anything special.
    ! they are named xtra1 .. xtra30, ixtra1 .. ixtra30, and lxtra1 .. lxtra30.
    ! they are automatically versioned, that is if you set s% xtra1, then
    ! s% xtra1_old will contains the value of s% xtra1 from the previous step
    ! and s% xtra1_older contains the one from two steps ago.

    if (s% log_center_density .lt. 9.0) return
    s% xtra1 = s% log_center_density

    ! this expression will evaluate to true if f times the log center density
    ! has crossed an integer during the last step.  If f = 5, then we will get
    ! output at log center density = {... 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ... }
    if ((floor(f * s% xtra1_old) - floor(f * s% xtra1) .ne. 0)) then

       ! save a profile & update the history
       s% need_to_update_history_now = .true.
       s% need_to_save_profiles_now = .true.

       ! by default the priority is 1; you can change that if you'd like
       s% save_profiles_model_priority = 3

    endif

    ! see extras_check_model for information about custom termination codes
    ! by default, indicate where (in the code) MESA terminated
    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

  end function extras_finish_step


  subroutine extras_after_evolve(id, id_extra, ierr)
    integer, intent(in) :: id, id_extra
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve


  ! routines for saving and restoring extra data so can do restarts

  subroutine alloc_extra_info(s)
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_alloc)
  end subroutine alloc_extra_info


  subroutine unpack_extra_info(s)
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_get)
  end subroutine unpack_extra_info


  subroutine store_extra_info(s)
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_put)
  end subroutine store_extra_info


  subroutine move_extra_info(s,op)
    type (star_info), pointer :: s
    integer, intent(in) :: op

    integer :: i, j, num_ints, num_dbls, ierr

    i = 0
    ! call move_int or move_flg
    num_ints = i

    i = 0
    ! call move_dbl

    num_dbls = i

    if (op /= extra_info_alloc) return
    if (num_ints == 0 .and. num_dbls == 0) return

    ierr = 0
    call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in star_alloc_extras'
       write(*,*) 'alloc_extras num_ints', num_ints
       write(*,*) 'alloc_extras num_dbls', num_dbls
       stop 1
    end if

  contains

    subroutine move_dbl(dbl)
      real(dp) :: dbl
      i = i+1
      select case (op)
      case (extra_info_get)
         dbl = s% extra_work(i)
      case (extra_info_put)
         s% extra_work(i) = dbl
      end select
    end subroutine move_dbl

    subroutine move_int(int)
      integer :: int
      i = i+1
      select case (op)
      case (extra_info_get)
         int = s% extra_iwork(i)
      case (extra_info_put)
         s% extra_iwork(i) = int
      end select
    end subroutine move_int

    subroutine move_flg(flg)
      logical :: flg
      i = i+1
      select case (op)
      case (extra_info_get)
         flg = (s% extra_iwork(i) /= 0)
      case (extra_info_put)
         if (flg) then
            s% extra_iwork(i) = 1
         else
            s% extra_iwork(i) = 0
         end if
      end select
    end subroutine move_flg

  end subroutine move_extra_info


end module run_star_extras

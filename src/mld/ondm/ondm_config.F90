#include "../../MLD_MACROS.INC"

module module_gin_in_hard 
  use module_kind_variables, only: kind_double 
  implicit none

  contains 

  subroutine gin_in_hard (tmp_lat, tmp_at, tmp_imcell, tmp_xc, tmp_itypc)
    !call gin_in_hard (lat, at, imcell)
    integer, intent(inout)  :: tmp_lat(3)
    real(kind_double), dimension(3, 3), intent(inout)    :: tmp_at
    integer, intent(inout)  :: tmp_imcell
    real(kind_double), allocatable, dimension(:,:), intent(inout)  :: tmp_xc
    integer, allocatable, dimension(:), intent(inout)  :: tmp_itypc 

    tmp_lat(1:3) = (/1,                1,              1/)
    tmp_at(1:3, 1) =(/11.421200000000001,      0.000000000000000,         0.000000000000000 /) 
    tmp_at(1:3, 2) =(/0.000000000000000,     11.421200000000001,         0.000000000000000 /)
    tmp_at(1:3, 3) =(/0.000000000000000,      0.000000000000000,        11.421200000000001/)
    tmp_imcell = 127
    if (allocated(tmp_xc)) deallocate(tmp_xc) ;  allocate (tmp_xc(tmp_imcell, 3))
    if (allocated (tmp_itypc)) deallocate(tmp_itypc) ; allocate (tmp_itypc(tmp_imcell))
    
    tmp_itypc(1:tmp_imcell) = 1

   tmp_xc( 1 ,1:3)=(/0.003590701,    0.003590701,    0.003590701 /)   
   tmp_xc( 2 ,1:3)=(/0.246562466,    0.003539316,    0.003539316 /)   
   tmp_xc( 3 ,1:3)=(/0.246555477,    0.246555477,    0.003527664 /)   
   tmp_xc( 4 ,1:3)=(/0.003539316,    0.246562466,    0.003539316 /)   
   tmp_xc( 5 ,1:3)=(/0.003539316,    0.003539316,    0.246562466 /)   
   tmp_xc( 6 ,1:3)=(/0.246555477,    0.003527664,    0.246555477 /)   
   tmp_xc( 7 ,1:3)=(/0.246536402,    0.246536402,    0.246536402 /)   
   tmp_xc( 8 ,1:3)=(/0.003527664,    0.246555477,    0.246555477 /)   
   tmp_xc( 9 ,1:3)=(/0.378458297,    0.125034628,    0.125034628 /)   
   tmp_xc(10 ,1:3)=(/0.374564305,    0.374564305,    0.125009370 /)   
   tmp_xc(11 ,1:3)=(/0.125034628,    0.378458297,    0.125034628 /)   
   tmp_xc(12 ,1:3)=(/0.125034628,    0.125034628,    0.378458297 /)   
   tmp_xc(13 ,1:3)=(/0.374564305,    0.125009370,    0.374564305 /)   
   tmp_xc(14 ,1:3)=(/0.373795381,    0.373795381,    0.373795381 /)   
   tmp_xc(15 ,1:3)=(/0.125009370,    0.374564305,    0.374564305 /)   
   tmp_xc(16 ,1:3)=(/0.000199014,    0.000199014,    0.500054449 /)   
   tmp_xc(17 ,1:3)=(/0.125038460,    0.125038460,    0.625035902 /)   
   tmp_xc(18 ,1:3)=(/0.000162622,    0.000162622,    0.750017235 /)   
   tmp_xc(19 ,1:3)=(/0.125070351,    0.125070351,    0.871599514 /)   
   tmp_xc(20 ,1:3)=(/0.000202509,    0.250230893,    0.500064301 /)   
   tmp_xc(21 ,1:3)=(/0.125006842,    0.375038706,    0.625018269 /)   
   tmp_xc(22 ,1:3)=(/0.000181590,    0.250251274,    0.749997716 /)   
   tmp_xc(23 ,1:3)=(/0.125027456,    0.374561548,    0.875472405 /)   
   tmp_xc(24 ,1:3)=(/0.000199014,    0.500054449,   -0.000199014 /)   
   tmp_xc(25 ,1:3)=(/0.125038460,    0.625035902,    0.125038460 /)   
   tmp_xc(26 ,1:3)=(/0.000202509,    0.500064301,    0.250230893 /)   
   tmp_xc(27 ,1:3)=(/0.125006842,    0.625018269,    0.375038706 /)   
   tmp_xc(28 ,1:3)=(/0.000022592,    0.499858795,    0.499858795 /)   
   tmp_xc(29 ,1:3)=(/0.124992316,    0.625034496,    0.625034496 /)   
   tmp_xc(30 ,1:3)=(/0.000008787,    0.499868390,    0.750191806 /)   
   tmp_xc(31 ,1:3)=(/0.125023232,    0.625032371,    0.875014367 /)   
   tmp_xc(32 ,1:3)=(/0.000162622,    0.750017235,   -0.000162622 /)   
   tmp_xc(33 ,1:3)=(/0.125070351,    0.871599514,    0.125070351 /)   
   tmp_xc(34 ,1:3)=(/0.000181590,    0.749997716,    0.250251274 /)   
   tmp_xc(35 ,1:3)=(/0.125027456,    0.875472405,    0.374561548 /)   
   tmp_xc(36 ,1:3)=(/0.000008787,    0.750191806,    0.499868390 /)   
   tmp_xc(37 ,1:3)=(/0.125023232,    0.875014367,    0.625032371 /)   
   tmp_xc(38 ,1:3)=(/0.000017821,    0.750206942,    0.750206942 /)   
   tmp_xc(39 ,1:3)=(/0.125051058,    0.875501314,    0.875501314 /)   
   tmp_xc(40  ,1:3)=(/0.250230893,    0.000202509,    0.500064301 /)   
   tmp_xc(41  ,1:3)=(/0.375038706,    0.125006842,    0.625018269 /)   
   tmp_xc(42  ,1:3)=(/0.250251274,    0.000181590,    0.749997716 /)   
   tmp_xc(43  ,1:3)=(/0.374561548,    0.125027456,    0.875472405 /)   
   tmp_xc(44  ,1:3)=(/0.250242878,    0.250242878,    0.500061816 /)   
   tmp_xc(45  ,1:3)=(/0.375041117,    0.375041117,    0.625024303 /)   
   tmp_xc(46  ,1:3)=(/0.250248664,    0.250248664,    0.749973942 /)   
   tmp_xc(47  ,1:3)=(/0.373800426,    0.373800426,    0.876242895 /)   
   tmp_xc(48  ,1:3)=(/0.250230893,    0.500064301,   -0.000202509 /)   
   tmp_xc(49  ,1:3)=(/0.375038706,    0.625018269,    0.125006842 /)   
   tmp_xc(50  ,1:3)=(/0.250242878,    0.500061816,    0.250242878 /)   
   tmp_xc(51  ,1:3)=(/0.375041117,    0.625024303,    0.375041117 /)   
   tmp_xc(52  ,1:3)=(/0.250044829,    0.499872297,    0.499872297 /)   
   tmp_xc(53  ,1:3)=(/0.375040693,    0.625041732,    0.625041732 /)   
   tmp_xc(54  ,1:3)=(/0.250023191,    0.499871416,    0.750166193 /)   
   tmp_xc(55  ,1:3)=(/0.375047373,    0.625025292,    0.874983671 /)   
   tmp_xc(56  ,1:3)=(/0.250251274,    0.749997716,   -0.000181590 /)   
   tmp_xc(57  ,1:3)=(/0.374561548,    0.875472405,    0.125027456 /)   
   tmp_xc(58  ,1:3)=(/0.250248664,    0.749973942,    0.250248664 /)   
   tmp_xc(59  ,1:3)=(/0.373800426,    0.876242895,    0.373800426 /)   
   tmp_xc(60  ,1:3)=(/0.250023191,    0.750166193,    0.499871416 /)   
   tmp_xc(61  ,1:3)=(/0.375047373,    0.874983671,    0.625025292 /)   
   tmp_xc(62  ,1:3)=(/0.250018476,    0.750175416,    0.750175416 /)   
   tmp_xc(63  ,1:3)=(/0.373794669,    0.876242569,    0.876242569 /)   
   tmp_xc(64  ,1:3)=(/0.500054449,    0.000199014,   -0.000199014 /)   
   tmp_xc(65  ,1:3)=(/0.625035902,    0.125038460,    0.125038460 /)   
   tmp_xc(66  ,1:3)=(/0.500064301,    0.000202509,    0.250230893 /)   
   tmp_xc(67  ,1:3)=(/0.625018269,    0.125006842,    0.375038706 /)   
   tmp_xc(68  ,1:3)=(/0.499858795,    0.000022592,    0.499858795 /)   
   tmp_xc(69  ,1:3)=(/0.625034496,    0.124992316,    0.625034496 /)   
   tmp_xc(70  ,1:3)=(/0.499868390,    0.000008787,    0.750191806 /)   
   tmp_xc(71  ,1:3)=(/0.625032371,    0.125023232,    0.875014367 /)   
   tmp_xc(72  ,1:3)=(/0.500064301,    0.250230893,   -0.000202509 /)   
   tmp_xc(73  ,1:3)=(/0.625018269,    0.375038706,    0.125006842 /)   
   tmp_xc(74  ,1:3)=(/0.500061816,    0.250242878,    0.250242878 /)   
   tmp_xc(75  ,1:3)=(/0.625024303,    0.375041117,    0.375041117 /)   
   tmp_xc(76  ,1:3)=(/0.499872297,    0.250044829,    0.499872297 /)   
   tmp_xc(77  ,1:3)=(/0.625041732,    0.375040693,    0.625041732 /)   
   tmp_xc(78  ,1:3)=(/0.499871416,    0.250023191,    0.750166193 /)   
   tmp_xc(79  ,1:3)=(/0.625025292,    0.375047373,    0.874983671 /)   
   tmp_xc(80  ,1:3)=(/0.499858795,    0.499858795,   -0.000022592 /)   
   tmp_xc(81  ,1:3)=(/0.625034496,    0.625034496,    0.124992316 /)   
   tmp_xc(82  ,1:3)=(/0.499872297,    0.499872297,    0.250044829 /)   
   tmp_xc(83  ,1:3)=(/0.625041732,    0.625041732,    0.375040693 /)   
   tmp_xc(84  ,1:3)=(/0.499678105,    0.499678105,    0.499678105 /)   
   tmp_xc(85  ,1:3)=(/0.625043802,    0.625043802,    0.625043802 /)   
   tmp_xc(86  ,1:3)=(/0.499672695,    0.499672695,    0.750356930 /)   
   tmp_xc(87  ,1:3)=(/0.625036538,    0.625036538,    0.874978068 /)   
   tmp_xc(88  ,1:3)=(/0.499868390,    0.750191806,   -0.000008787 /)   
   tmp_xc(89  ,1:3)=(/0.625032371,    0.875014367,    0.125023232 /)   
   tmp_xc(90  ,1:3)=(/0.499871416,    0.750166193,    0.250023191 /)   
   tmp_xc(91  ,1:3)=(/0.625025292,    0.874983671,    0.375047373 /)   
   tmp_xc(92  ,1:3)=(/0.499672695,    0.750356930,    0.499672695 /)   
   tmp_xc(93  ,1:3)=(/0.625036538,    0.874978068,    0.625036538 /)   
   tmp_xc(94  ,1:3)=(/0.499687124,    0.750373882,    0.750373882 /)   
   tmp_xc(95  ,1:3)=(/0.625042998,    0.875000413,    0.875000413 /)   
   tmp_xc(96  ,1:3)=(/0.750017235,    0.000162622,   -0.000162622 /)   
   tmp_xc(97  ,1:3)=(/0.871599514,    0.125070351,    0.125070351 /)   
   tmp_xc(98  ,1:3)=(/0.749997716,    0.000181590,    0.250251274 /)   
   tmp_xc(99  ,1:3)=(/0.875472405,    0.125027456,    0.374561548 /)   
   tmp_xc(100 ,1:3)=(/0.750191806,    0.000008787,    0.499868390 /)   
   tmp_xc(101 ,1:3)=(/0.875014367,    0.125023232,    0.625032371 /)   
   tmp_xc(102 ,1:3)=(/0.750206942,    0.000017821,    0.750206942 /)   
   tmp_xc(103 ,1:3)=(/0.875501314,    0.125051058,    0.875501314 /)   
   tmp_xc(104 ,1:3)=(/0.749997716,    0.250251274,   -0.000181590 /)   
   tmp_xc(105 ,1:3)=(/0.875472405,    0.374561548,    0.125027456 /)   
   tmp_xc(106 ,1:3)=(/0.749973942,    0.250248664,    0.250248664 /)   
   tmp_xc(107 ,1:3)=(/0.876242895,    0.373800426,    0.373800426 /)   
   tmp_xc(108 ,1:3)=(/0.750166193,    0.250023191,    0.499871416 /)   
   tmp_xc(109 ,1:3)=(/0.874983671,    0.375047373,    0.625025292 /)   
   tmp_xc(110 ,1:3)=(/0.750175416,    0.250018476,    0.750175416 /)   
   tmp_xc(111 ,1:3)=(/0.876242569,    0.373794669,    0.876242569 /)   
   tmp_xc(112 ,1:3)=(/0.750191806,    0.499868390,   -0.000008787 /)   
   tmp_xc(113 ,1:3)=(/0.875014367,    0.625032371,    0.125023232 /)   
   tmp_xc(114 ,1:3)=(/0.750166193,    0.499871416,    0.250023191 /)   
   tmp_xc(115 ,1:3)=(/0.874983671,    0.625025292,    0.375047373 /)   
   tmp_xc(116 ,1:3)=(/0.750356930,    0.499672695,    0.499672695 /)   
   tmp_xc(117 ,1:3)=(/0.874978068,    0.625036538,    0.625036538 /)   
   tmp_xc(118 ,1:3)=(/0.750373882,    0.499687124,    0.750373882 /)   
   tmp_xc(119 ,1:3)=(/0.875000413,    0.625042998,    0.875000413 /)   
   tmp_xc(120,1:3)=(/0.750206942,    0.750206942,    0.000017821 /)   
   tmp_xc(121,1:3)=(/0.875501314,    0.875501314,    0.125051058 /)   
   tmp_xc(122,1:3)=(/0.750175416,    0.750175416,    0.250018476 /)   
   tmp_xc(123,1:3)=(/0.876242569,    0.876242569,    0.373794669 /)   
   tmp_xc(124,1:3)=(/0.750373882,    0.750373882,    0.499687124 /)   
   tmp_xc(125,1:3)=(/0.875000413,    0.875000413,    0.625042998 /)   
   tmp_xc(126,1:3)=(/0.750389797,    0.750389797,    0.750389797 /)   
   tmp_xc(127,1:3)=(/0.876289622,    0.876289622,    0.876289622 /)  
  end subroutine gin_in_hard

  subroutine gin_in_soft_get_imm 
    use, intrinsic :: iso_fortran_env, dp=>real64
    use ondm_gen_com_m, only: imm, im, imm_glob, im_glob, fnam  
    use mld_logger
    implicit none 
    logical :: ginok 
    integer :: lugin 
    character(len=:), allocatable    :: fnamgin
    integer :: la, lb, lc, imcell  
    real(dp) :: at(3,3) 
    fnamgin = fnam//'.gin'
    inquire (file=fnamgin, exist=ginok)
    if (ginok) then 
      open (newunit=lugin, file=fnamgin, status='unknown')
      call log_info('NDM -- BUILD LATTICE  FROM  -- '//trim(fnamgin))
    else 
      return   
    end if 
    
    if (ginok) then 
       read (lugin, *) la, lb, lc
       ! a
       read (lugin, *) at(1, 1), at(2, 1), at(3, 1)
       ! b
       read (lugin, *) at(1, 2), at(2, 2), at(3, 2)
       ! c
       read (lugin, *) at(1, 3), at(2, 3), at(3, 3) 
       read (lugin, *) imcell  ! number of atoms in UC

       if (imcell > imm_glob) then 
         imm_glob = imcell 
         imm = imcell 
         im_glob = imcell
         im = imcell 
       end if
       im_glob = la*lb*lc*imcell

       if (im_glob > imcell) then 
         imm_glob = im_glob
         im = im_glob
         imm = im_glob
       end if

       close(lugin, status='keep')

    end if   !if ginok  


  end subroutine gin_in_soft_get_imm 
end module module_gin_in_hard 

subroutine ondm_config()
  !********************************************
  !    read gin configuration file at minimum. 
  !********************************************

  USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m, only: imm, im, imm_glob, im_glob, fnam, lat, at, lvpread, & 
                       normat, zl, zls2, bg, lperiod, nzl, lsuivinonpbc, A2cm
  use ondm_var_pot, only: rumax, alpha, na 
  use ondm_tab_imm_m, only: xp, xpnonpbc, ityp, num_at_glob, ax, axnonpbc 
  use module_gin_in_hard, only: gin_in_hard
  use ondm_transform_coord, only: ondm_cryst_to_cart, ondm_period, ondm_recips


  use mld_logger
  use mld_string
  use mld_mpi

  implicit none

  character(len=:), allocatable    :: fnamgin
  integer  :: i, ia, ib, ic, icell, &
              lugin, imcell, la, lb, lc

  integer, dimension(:), pointer   :: itypc
  real(kind_double), dimension(:, :), pointer :: xc
  integer, dimension(:), pointer   :: ibuffer
  real(kind_double), dimension(:, :), pointer :: buffer

  integer  :: tmp_lat(3)
  real(kind_double), dimension(3, 3)    :: tmp_at
  integer :: tmp_imcell
  real(kind_double), allocatable, dimension(:,:) :: tmp_xc
  integer, allocatable, dimension(:)  :: tmp_itypc 

  real(kind_double)   :: rumax_init
  real(kind_double)   :: alpha_init
  ! integer , dimension(imm,ntyp) :: fv
  ! integer , dimension(6000,10) :: fv    ! Truc_bizarre_jmd
  logical :: ginok

  !-----------------------------------------------------
  ! READING FROM THE CONFIGURATION FILE
  !---------------------------------------------------
  allocate (ibuffer(imm_glob))
  allocate (buffer(3, imm_glob))


  lvpread = .false.
  ! open fichier .gin
  !lugin = 92
  fnamgin = fnam//'.gin'
  inquire (file=fnamgin, exist=ginok)
  if (ginok) then 
    open (newunit=lugin, file=fnamgin, status='unknown')
    call log_info('ML: oNDM -- BUILD LATTICE  FROM  -- '//trim(fnamgin))
  else 
    call gin_in_hard(tmp_lat, tmp_at, tmp_imcell, tmp_xc, tmp_itypc)
    !debug call log_info('NDM -- BUILD LATTICE  FROM  HARD GIN ! ')
  end if   


  ! number of cells in 3 directions
  if (ginok) then 
    read (lugin, *) lat(1), lat(2), lat(3)
  else 
    lat(1:3) = tmp_lat(1:3)
  end if

  if (ginok) then 
    ! **** coordonnes des vecteurs de maille en A dans une base orthonormee ****
    ! a
    read (lugin, *) at(1, 1), at(2, 1), at(3, 1)
    ! b
    read (lugin, *) at(1, 2), at(2, 2), at(3, 2)
    ! c
    read (lugin, *) at(1, 3), at(2, 3), at(3, 3)
  else 
    at(1:3, 1:3) = tmp_at(1:3, 1:3)
  end if 
  
  do ic = 1, 3
    normat(ic) = 0
    !TODOcmA units 
    at(:, ic) = at(:, ic) * A2cm * lat(ic)
    normat(ic) = normat(ic) + sum(at(:, ic)**2)
    normat(ic) = sqrt(normat(ic))
    zl(ic) = normat(ic)
    zls2(ic) = zl(ic)/2.
  end do
  la = lat(1)
  lb = lat(2)
  lc = lat(3)

  ! Il est important de conserver rue et alpha identique a
  ! chaque appel a la routine divid, on sauvegarde donc la valeur
  ! initiale pour la remettre en sortie
  rumax_init = rumax
  alpha_init = alpha
  call ondm_recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
  do ic = 1, 3
    normat(ic) = sqrt(sum(bg(:, ic)**2))
    nzl(ic) = 1.0/normat(ic)
    ! if(rang==0)  write (6,)'nzl',nzl(ic)*1d8
  end do

  
  call ondm_divid(0)
  rumax = rumax_init
  alpha = alpha_init

  !debug  call log_debug('NDM -- read imcell')
  if (ginok) then 
    read (lugin, *) imcell  ! number of atoms in UC
  else 
    imcell = tmp_imcell 
  end if 
  !debug  call log_info('NDM -- atoms per unit cell '//vtoa(imcell))
  if (imcell > imm_glob) call mld_mpi_abort('too much atoms in unit cell ()')

  allocate (xc(imcell, 3))
  allocate (itypc(imcell))


  im_glob = la*lb*lc*imcell
  if (im_glob > imm_glob) call mld_mpi_abort('imm_glob too small')
  !debug call log_info('NDM -- created cells number '//vtoa([la, lb, lc, la*lb*lc]))

  do i = 1, imcell
    ! write (6,*)i,imcell
    if (ginok) then 
    read (lugin, *) xc(i, 1), xc(i, 2), xc(i, 3), itypc(i)
    else 
      itypc(i) = tmp_itypc(i) 
      xc(i,1:3) = tmp_xc(i,1:3)
    end if   
  end do


  if (lperiod .EQV. .true.) then
    do i = 1, imcell
      WHERE ((xc(i, :) .LT. 0.d0) .OR. (xc(i, :) .GE. 1.d0))
        xc(i, :) = xc(i, :) - Dble(Floor(xc(i, :)))
      END WHERE
    end do
  end if

  i = 0
  im = 0


  do ia = 1, la
    do ib = 1, lb
      do ic = 1, lc
        do icell = 1, imcell
          i = i + 1
          im = im + 1
          xp(1, i) = (xc(icell, 1) + float(ia - 1))/float(la)
          xp(2, i) = (xc(icell, 2) + float(ib - 1))/float(lb)
          xp(3, i) = (xc(icell, 3) + float(ic - 1))/float(lc)

          ityp(i) = itypc(icell)
          num_at_glob(i) = i
        end do
      end do
    end do
  end do

  deallocate (xc)
  deallocate (itypc)

  na = 0
  do i = 1, im
    na(ityp(i)) = na(ityp(i)) + 1
    ax(:, i) = xp(:, i)
  end do

  call ondm_cryst_to_cart(imm, xp, at, 1)               ! cryst vers cart
  ax(:, :im) = xp(:, :im)
  if ((lperiod) .and. (lsuivinonpbc)) then
    call ondm_cryst_to_cart(imm, xpnonpbc, at, 1)         ! cryst vers cart
    axnonpbc(:, :im) = xpnonpbc(:, :im)
  end if


  if (ginok) close (lugin)

  if (lperiod .EQV. .true.) call ondm_period

  !$! if (rang == 0) then
  !$!   tmp = '-- SIMULATION BOX --'//nwl// &
  !$!         'atoms per box '//vtoa(im_glob)//nwl// &
  !$!         'box size ZL '//vtoa([1D+08*zl(1), 1D+08*zl(2), 1D+08*zl(3)])//nwl
  !$!   do i = 1, 3
  !$!     tmp = tmp//'vector    '//vtoa(i)//vtoa(at(:, i)*1.0d8)//nwl
  !$!   end do
  !$!   do iti = 1, ntyp
  !$!     if (na(iti) == 0) cycle
  !$!     tmp = tmp//vtoa(na(iti))//'atoms of type'//vtoa(iti)//nwl
  !$!   end do
  !$!   call log_debug(tmp)
  !$! end if

  !$! if (llangevin .eqv. .true.) allocate (Gl(3, imm))
  deallocate (ibuffer)
  deallocate (buffer)

  return


end subroutine ondm_config




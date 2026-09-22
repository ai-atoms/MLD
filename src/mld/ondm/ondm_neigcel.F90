
#include "../../MLD_MACROS.INC"

subroutine ondm_neigcel
  ! subroutine cells forwards cells
  ! version du 28 septembre 2000

  !USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m, only: nox, noy, noz, noxyz, ncel, deltadist
  use mld_logger

  implicit none
  integer  :: kx, ky, kz, koo, l, lz, mz, ly, my, lx, mx, kxy

  _NAMECURRENT_("neigcel")
  _MLD_BEGIN_ 


  if (noxyz == 1) then
    ncel(1, 0) = 1
    deltadist = 0
  else

    do kz = 1, noz
      do ky = 1, noy
        do kx = 1, nox
          koo = 1 + (kx - 1) + nox*((ky - 1) + noy*(kz - 1))

          ncel(koo, 0) = koo
          deltadist(:, 0, koo) = 0

          l = 1
          do lz = -1, 1
            do ly = -1, 1
              do lx = -1, 1
                deltadist(:, l, koo) = 0

                mz = kz + lz
                if (mz < 1) then
                  mz = mz + noz
                  deltadist(3, l, koo) = 1
                end if
                if (mz > noz) then
                  mz = mz - noz
                  deltadist(3, l, koo) = -1
                end if

                my = ky + ly
                if (my < 1) then
                  my = my + noy
                  deltadist(2, l, koo) = 1
                end if
                if (my > noy) then
                  my = my - noy
                  deltadist(2, l, koo) = -1
                end if

                mx = kx - lx
                if (mx < 1) then
                  mx = mx + nox
                  deltadist(1, l, koo) = 1
                end if
                if (mx > nox) then
                  mx = mx - nox
                  deltadist(1, l, koo) = -1
                end if

                kxy = 1 + (mx - 1) + nox*((my - 1) + noy*(mz - 1))
                if (kxy == koo) cycle
                ncel(koo, l) = kxy
                !                        write (6,*)koo,lz,ly,lx,l,kxy
                !                        if ((kz==noz).and.(lz==1))write (6,*)koo,lz,l,kxy
                !                        if ((kz==1).and.(lz==-1))write (6,*)koo,lz,l,kxy
                l = l + 1
              end do
            end do
          end do
        end do
      end do
    end do
  end if
  _MLD_END_

end subroutine ondm_neigcel

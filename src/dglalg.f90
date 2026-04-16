subroutine dglalg_f(x_l, x_r, z_l, z_r, n, maxtime, ctime, &
                    x_val, nx, t_val, nt,                    &
                    w_out, f_out, counter, info)

  implicit none

  ! --- inputs ---
  integer,          intent(in)  :: n, nx, nt
  double precision, intent(in)  :: x_l(n), x_r(n), z_l(n), z_r(n)
  double precision, intent(in)  :: maxtime, ctime
  double precision, intent(in)  :: x_val(nx), t_val(nt)

  ! --- outputs ---
  double precision, intent(out) :: w_out(nx), f_out(nt)
  integer,          intent(out) :: counter, info

  ! --- local ---
  integer          :: i, j, k, iter
  double precision :: w_new(nx), w_old(nx)
  double precision :: f_new(nt), f_old(nt)
  double precision :: w_acc(nx), f_acc(nt)
  double precision :: denom, contrib, tol, diff
  logical          :: hit
  double precision, parameter :: eps = 0.0005d0

  info    = 0
  counter = 0

  ! initialize
  do j = 1, nx
    w_new(j) = 1.0d0 / dble(nx)
  end do
  do k = 1, nt
    f_new(k) = 1.0d0 / dble(nt)
  end do

  tol = 1.0d10

  do while (tol > eps)

    if (counter == 1000) then
      info = 1
      return
    end if

    w_old = w_new
    f_old = f_new

    ! zero accumulators
    do j = 1, nx
      w_acc(j) = 0.0d0
    end do
    do k = 1, nt
      f_acc(k) = 0.0d0
    end do

    ! main pass over subjects
    do i = 1, n

      ! compute denominator for subject i
      denom = 0.0d0
      do k = 1, nt
        do j = 1, nx
          if (x_l(i) <= x_val(j) .and. x_val(j) <= x_r(i) .and. &
              z_l(i) <= x_val(j) + t_val(k) .and.                &
              x_val(j) + t_val(k) <= z_r(i)) then
            denom = denom + w_new(j) * f_new(k)
          end if
        end do
      end do

      if (denom == 0.0d0) then
        info = -1
        return
      end if

      ! accumulate w and f updates
      do k = 1, nt
        do j = 1, nx
          if (x_l(i) <= x_val(j) .and. x_val(j) <= x_r(i) .and. &
              z_l(i) <= x_val(j) + t_val(k) .and.                &
              x_val(j) + t_val(k) <= z_r(i)) then
            contrib      = w_new(j) * f_new(k) / denom
            w_acc(j)     = w_acc(j) + contrib
            f_acc(k)     = f_acc(k) + contrib
          end if
        end do
      end do

    end do  ! end subject loop

    ! normalize
    do j = 1, nx
      w_new(j) = w_acc(j) / dble(n)
    end do
    do k = 1, nt
      f_new(k) = f_acc(k) / dble(n)
    end do

    ! compute tolerance
    tol = 0.0d0
    do j = 1, nx
      diff = w_old(j) - w_new(j)
      if (diff < 0.0d0) diff = -diff
      tol = tol + diff
    end do
    do k = 1, nt
      diff = f_old(k) - f_new(k)
      if (diff < 0.0d0) diff = -diff
      tol = tol + diff
    end do

    counter = counter + 1

  end do  ! end EM loop

  w_out = w_new
  f_out = f_new

end subroutine dglalg_f
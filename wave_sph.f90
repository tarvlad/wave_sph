module constants
implicit none
real, parameter :: pi = acos(-1.0)
real, parameter  :: rho_0 = 1.0
real, parameter  :: delta = 1.0 / (10.0**3)
real, parameter  :: cs = 1.0
real, parameter :: l_bound = -0.7
real, parameter  :: r_bound = 1.7
integer, parameter  :: n_particles_per_unit_side = 100
real, parameter  :: h = 0.1
real, parameter  :: time_moment = 1.0
real, parameter  :: tau = 0.00001
real, parameter  :: k = 2 * pi
real, parameter  :: m_g = rho_0 / n_particles_per_unit_side
integer, parameter :: n_g = n_particles_per_unit_side * (r_bound - l_bound)

end module constants


real function rho_analytical(x, t)
use constants
implicit none
real, intent(in) :: x
real, intent(in) :: t
rho_analytical = rho_0 + delta * cos(k * x - k * cs * t)
end function rho_analytical


real function v_analytical(x, t)
use constants
implicit none
real, intent(in) :: x
real, intent(in) :: t
v_analytical = delta * (cs / rho_0) * cos(k * x - k * cs * t)
end function v_analytical


subroutine init_particles(
  p_coords, p_velocities, p_densities,
  r_p_coords, r_p_velocities, r_p_densities
)
use constants
implicit none
real, dimension(n_g), intent(inout) :: p_coords
real, dimension(n_g), intent(inout) :: p_velocities
real, dimension(n_g), intent(inout) :: p_densities
real, dimension(n_g), intent(inout) :: r_p_coords
real, dimension(n_g), intent(inout) :: r_p_velocities
real, dimension(n_g), intent(inout) :: r_p_densities

real, dimension(n_g + 1) :: borders
real :: delta_b
real :: coord
real :: curr_b
real :: curr_i
real :: velocity
integer :: i

borders(1) = l_bound
borders(n_g + 1) = r_bound

delta_b = (r_bound - l_bound) / ((10.0**5) * real(n_g))
coord = 0.0

do i = 1, n_g - 1
  curr_b = borders(i)
  curr_i = 0.0

  do while (curr_i < m_g)
    curr_b = curr_b + delta_b
    coord = curr_b
    curr_i = curr_i + delta_b * rho_analytical(coord, 0.0)
  end do

  borders(i + 1) = curr_b
end do

do i = 1, n_g
  coord = (borders(i) + borders(i + 1)) / 2.0
  p_coords(i) = coord
end do

do i = 1, n_g
  velocity = v_analytical(p_coords(i), 0.0)
  p_velocities(i) = velocity
  p_densities(i) = rho_analytical(p_coords(i), 0.0)

  r_p_coords(i) = 0.0
  r_p_densities(i) = 0.0
  r_p_velocities(i) = 0.0
end do

end subroutine init_particles


real function init_rho_max_error(p_coords, p_densities)
use constants
implicit none
real, dimension(n_g), intent(inout) :: p_coords
real, dimension(n_g), intent(inout) :: p_densities

integer :: i
real :: error
real :: max_error
error = 0.0
max_error = 0.0

do i = 1, n_g
  if ((p_coords(i).ge.0.0).and.(p_coords(i).le.1.0)) then
    error = abs(p_densities(i) - rho_analytical(p_coords(i), 0.0))
    if (error.gt.max_error) then
      max_error = error
    end if
  end if 
end do

init_rho_max_error = max_error
end function init_rho_max_error


real function w(x_a, x_b)
use constants
implicit none
real, intent(in) :: x_a
real, intent(in) :: x_b

real :: q
real :: y
q = abs(x_a - x_b) / h

if (q.le.1.0) then !S4O2
  y = 3.0 / (2.0 * h)
  y = y * (1.0 - q)**5
  y = y * (8 * q**2 + 5.0 * q + 1.0)
else
  y = 0.0
end if

w = y
end function w


real function dw(x_a, x_b)
use constants
implicit none
real, intent(in) :: x_a
real, intent(in) :: x_b

real :: q
real :: y
q = abs(x_a - x_b) / h

if (q.le.1.0) then !S4O2
  y = 3.0 / (2.0 * h)
  y = y * (1.0 - q)**4
  y = y * (-56.0 * q**2 - 14 * q)
else
  y = 0.0
end if

dw = y
end function dw


real function w_grad(x_a, x_b)
use constants
implicit none
real, intent(in) :: x_a
real, intent(in) :: x_b

real :: r
r = x_a - x_b

if (abs(r).ne.0.0) then
  w_grad = r / (abs(r) * h) * dw(x_a, x_b)
else 
  w_grad = 0.0
end if

end function w_grad


subroutine upd(
  p_coords, p_velocities, p_densities,
  r_p_coords, r_p_velocities, r_p_densities
)
use constants
implicit none
real, dimension(n_g), intent(inout) :: p_coords
real, dimension(n_g), intent(inout) :: p_velocities
real, dimension(n_g), intent(inout) :: p_densities

real, dimension(n_g), intent(inout) :: r_p_coords
real, dimension(n_g), intent(inout) :: r_p_velocities
real, dimension(n_g), intent(inout) :: r_p_densities

real :: r_ij
integer :: i
integer :: j

do i = 1, n_g
  r_p_coords(i) = r_p_velocities(i)
  r_p_coords(i) = p_coords(i) + tau * r_p_coords(i)

  r_p_densities(i) = 0.0
  do j = 1, n_g
    r_ij = p_coords(i) - p_coords(j)
    if (abs(r_ij).le.h) then
      r_p_densities(i) = r_p_densities(i) + &
        (p_velocities(i) + p_velocities(j)) * w_grad(p_coords(i), p_coords(j))
    end if
  end do
  r_p_densities(i) = r_p_densities(i) * (-m_g)
  r_p_densities(i) = r_p_densities(i) * tau + p_densities(i)

  r_p_velocities(i) = 0.0
  do j = 1, n_g
    r_ij = p_coords(i) - p_coords(j)
    if (abs(r_ij).le.h) then
      r_p_velocities(i) = r_p_velocities(i) + w_grad(p_coords(i), p_coords(j))
    end if
  end do
  r_p_velocities(i) = r_p_velocities(i) * (-1.0) * cs**2 * m_g / p_densities(i)
  r_p_velocities(i) = p_velocities(i) + tau * r_p_velocities(i)
end do

do i = 1, n_g
  p_coords(i) = r_p_coords(i)
  p_velocities(i) = r_p_velocities(i)
  p_densities(i) = r_p_densities(i)
end do

end subroutine upd


program wave_sph
use constants
implicit none
real, dimension(n_g) :: p_coords
real, dimension(n_g) :: p_velocities
real, dimension(n_g) :: p_densities

real, dimension(n_g) :: r_p_coords
real, dimension(n_g) :: r_p_velocities
real, dimension(n_g) :: r_p_densities

real :: init_max_error_rho
integer :: n_time
integer :: i
real :: max_err_rho
real :: max_err_v
real :: err_tmp

call init_particles(
  p_coords, p_velocities, p_densities,
  r_p_coords, r_p_velocities, r_p_densities
)

init_max_error_rho = init_rho_max_error(p_coords, p_densities)
n_time = ceiling(time_moment / tau)

do i = 1, n_g
  call upd(
    p_coords, p_velocities, p_densities,
    r_p_coords, r_p_velocities, r_p_densities
  )
end do

max_err_rho = 0.0
max_err_v = 0.0
do i = 1, n_g
  if ((p_coords(i).ge.0.0).and.(p_coords(i).le.1.0)) then
    err_tmp = abs(p_velocities(i) - v_analytical(p_coords(i), time_moment))
    if (err_tmp.gt.max_err_v) then
      max_err_v = err_tmp
    end if
    err_tmp = abs(p_densities(i) - rho_analytical(p_coords(i), time_moment))
    if (err_tmp.gt.max_err_rho) then
      max_err_rho = err_tmp
    end if
  end if
end do

write(*, *) "tau ", tau
write(*,*) "h ", h
write(*,*) "N ", n_particles_per_unit_side
write(*,*) "time ", time_moment
write(*,*) ""
write(*,*) "max_v_err", max_err_v
write(*,*) "max_rho_err", max_err_rho
write(*,*) ""

write(*,*) "x"
write(*,*) "  v rho" 
write(*,*) "  a_v a_rho"
do i = 1, n_g
  if ((p_coords(i).ge.0.0).and.(p_coords(i).le.1.0)) then
    write(*,*) p_coords(i)
    write(*,*) "  ", p_velocities(i), " ", p_densities(i)
    write(*,*) "  ", v_analytical(p_coords(i), time_moment), " ", rho_analytical(p_coords(i), time_moment)
  end if
end do


end program wave_sph
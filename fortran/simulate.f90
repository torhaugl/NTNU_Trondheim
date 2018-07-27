module input
   implicit none

   ! Input file / Parameters
   real,    parameter :: dt = 0.07/60.0 !time step (min)
   real,    parameter :: final_time = 0.5*60.0 !minutes
   integer, parameter :: NumTrials = final_time / dt
   real,    parameter :: c_bulk = 0.2 ! Concentration bulk substrate
   real,    parameter :: v_width = 17.0 !Voxel length
   real,    parameter :: diff_s = 40680.0, diff_q = 33300.0 !Diffusion subst./QSM
   integer, parameter :: v_size(3) = (/10,10,100/)
   integer, parameter :: v_count = v_size(1) * v_size(2) * v_size(3)
   real,    parameter :: Vmax = 46e-3 !Maximum substrate uptake
   real,    parameter :: Ks = 2.34e-3 ! Half-saturaton const (substrate uptake)
   real,    parameter :: Zqd = 8.3, Zqu = 1230.0 !QSM production 
   real,    parameter :: Kq = 10
   real,    parameter :: Ymax = 0.444 
   real,    parameter :: maintenance = 0.6e-3 
   real,    parameter :: avg_mass_cell = 420.0
   real,    parameter :: density_cell = 290.0
   real,    parameter :: max_vol_frac = 0.52
   real,    parameter :: Mmax = 14700.0
   integer, parameter :: Nmax = floor(density_cell * max_vol_frac * v_width**3 / Mmax) ! Standard input: 50
   real,    parameter :: alpha = 1.33
   real,    parameter :: beta = 10
   real,    parameter :: gamma = 0.1
   real,    parameter :: eps_mass = avg_mass_cell
   real,    parameter :: Zed = 0, Zeu = 1e-3
   real,    parameter :: mu = 1e-3 !Transfer coefficient
end

program simulate
   ! Simulates time-steps
   use input
   implicit none

   ! Functions
   real, external :: avg

   ! Variables
   real,    dimension(v_count)         :: prod_s, prod_q
   integer, dimension(v_count,2)       :: eps_count
   real,    dimension(v_count,2)       :: c_s, c_q, up, eps_amount, pressure
   real,    dimension(Nmax,v_count,2)  :: biomass
   integer :: pos(3) !xyz position

   ! Loops
   integer :: i, n ! Standard indexes for loop
   integer :: i_list_neigh(6) !Index for neighbours of a voxel
   real :: start, finish

   ! Initialize arrays
   c_s(:,:) = c_bulk
   c_q(:,:) = 0.0
   up(:,:) = 0.0
   biomass(:,:,:) = 0.0
   eps_count(:,:) = 0
   eps_amount(:,:) = 0.0
   pressure(:,:) = 0


   do i = 1, v_size(1)*v_size(2)
      call append_biomass(600.0,i, biomass)
      biomass(:,i,1) = biomass(:,i,2)
   end do

   !!! CODE
   print*,"CODE START"
   print*,"Number of steps:", NumTrials
   print*,

   call cpu_time(start)

   ! TODO Updates independent of neighbours (all except concentration)
   !      can be calculated outside voxel loop
   do n = 1,NumTrials
      do i = 1, v_count
         call update_eps(i, biomass, up, eps_amount, eps_count)
         call update_mass(i, c_s, biomass)
         call update_stochastics(i,c_q,biomass,up)
         call update_pressure(i, biomass, eps_count, pressure)
         call update_displacement(i, pressure, biomass, eps_amount, eps_count)
         call update_production_s(i, c_s, biomass, prod_s)
         call update_production_q(i, c_q, biomass, up, prod_q)
         ! Remove above out of voxel loop

         call update_concentration(i, prod_s, diff_s, c_s) !Substrate
         call update_concentration(i, prod_q, diff_q, c_q) !QSM
      end do
      c_s(:,1) = c_s(:,2) ! Insert the newly calculated concentrations
      c_q(:,1) = c_q(:,2)
      biomass(:,:,1) = biomass(:,:,2)
      eps_count(:,1) = eps_count(:,2)
      eps_amount(:,1) = eps_amount(:,2)
      up(:,1) = up(:,2)
   end do

   call cpu_time(finish)

   print*,"CPU time(s):", finish-start
   print*,"Model time(min):", final_time
   print*,
   print*,"Top"   ,c_s(10000,1)  ,c_q(10000,1)  ,eps_amount(10000,1)
   print*,"Mid"   ,c_s(5000,1)   ,c_q(5000,1)   ,eps_amount(5000,1)
   print*,"Bot"   ,c_s(1,1)      ,c_q(1,1)      ,eps_amount(1,1)
   print*,"AVG"   ,avg(c_s(:,1)) ,avg(c_q(:,1)) ,avg(eps_amount(:,1))
   print*,"Bio"   ,biomass(:4,1,1)
   print*,"upp"   ,up(:4,1)
   call get_count_up(1,biomass,up,i)
   print*,"Up:"   ,i
   call get_count_down(1,biomass,up,i)
   print*,"Dow"   ,i
   call get_count_particles(biomass(:,1,1),eps_count(1,1),i)
   print*,"Par"   ,i

end program simulate



subroutine update_pressure(i, biomass, eps_count, pressure)
   use input
   implicit none

   integer, intent(in) :: i, eps_count(v_count,2)
   real, intent(in) :: biomass(Nmax,v_count,2)
   real, intent(out) :: pressure(v_count,2)
   integer :: count

   call get_count_particles(biomass(:,i,1), eps_count(i,1), count)
   if (count == Nmax) then
      pressure(i,2) = 1.0
   else 
      pressure(i,2) = count / (Nmax - count)
   end if

   !if (pressure(i,2) > 1) print*, "Pressure > 1"
end

subroutine update_displacement(i, pressure, biomass, eps_amount, eps_count)
   use input
   implicit none

   integer, intent(in) :: i
   real, intent(in) :: pressure(v_count,2)
   integer, intent(out) :: eps_count(v_count,2)
   real, intent(out) :: biomass(Nmax,v_count,2), eps_amount(v_count,2)

   integer :: j,k, j_list_neigh(6), particle_transfer, count, count_neigh


   call get_count_particles(biomass(:,i,1), eps_count(i,1), count)
   call get_index_neighbours(i, j_list_neigh)

   particle_transfer = 0
   do j =1,6
      k = j_list_neigh(j)
      if (.NOT. k == 0) then
         call get_count_particles(biomass(:,k,1), eps_count(k,1), count_neigh)
         particle_transfer = particle_transfer + floor(mu * (pressure(i,1) - pressure(k,1)) *(count - count_neigh) )
      end if
   end do

   if (particle_transfer < 0) print*, "Error, particle_transfer<0", particle_transfer
   if (particle_transfer > 0) print*, "Particle transfer", particle_transfer



end

real function avg(arr)
   ! Must be length v_count
   use input
   implicit none
   real, intent(in) :: arr(v_count)
   avg = sum(arr) / size(arr)
   return
end

subroutine update_eps(i, biomass, up, eps_amount, eps_count)
   ! Calculate new eps. Create particle if eps > eps_mass
   ! TODO Check total particle count (Pressure?)
   use input !v_count Nmax Zed Zeu eps_mass
   implicit none

   integer, intent(in) :: i
   real, dimension(Nmax, v_count,2), intent(in)  :: biomass
   real, dimension(v_count,2), intent(out) :: eps_amount, up
   integer, intent(out) :: eps_count(v_count,2)
   integer :: count_up, count_down

   call get_count_up(i,biomass,up,count_up)
   call get_count_down(i,biomass,up,count_down)

   eps_amount(i,2) = eps_amount(i,1) + Zed*count_down + Zeu*count_up

   if (eps_amount(i,2) >= eps_mass) then
      eps_amount(i,2) = eps_amount(i,2) - eps_mass
      eps_count(i,2) = eps_count(i,1) + 1

      if (eps_count(i,2) > Nmax) then
         print*, "EPS count too high", i, eps_count(i,2)
      end if
   end if
end

subroutine update_mass(i, c_s, biomass)
   use input !v_count Nmax Ymax, Ks, Vmax, m, dt
   implicit none

   integer, intent(in) :: i
   real, intent(in) :: c_s(v_count, 2)
   real, intent(out) :: biomass(Nmax,v_count,2)

   biomass(:,i,2) = biomass(:,i,1) &
      + dt*Ymax* ( c_s(i,2)/(Ks + c_s(i,2) ) - maintenance)*biomass(:,i,1)
end

subroutine update_production_q(i, c_q, biomass, up, prod_q)
   ! Output the production of QSM in voxel i
   use input !v_count,Nmax, Kq, Zqd, Zqu
   implicit none

   integer, intent(in)  :: i
   real, intent(in)     :: up(v_count,2), c_q(v_count,2)
   real, intent(in), dimension(Nmax,v_count,2) :: biomass
   real, intent(out)    :: prod_q(v_count)
   integer :: count_up, count_down

   call get_count_up(i,biomass,up,count_up)
   call get_count_down(i,biomass,up,count_down)

   prod_q(i) = Zqd*count_down + Zqu*count_up * c_q(i,1)/(Kq + c_q(i,1) )
end

subroutine get_count_up(i,biomass,up,count_up)
   ! TODO Worth it?
   use input
   implicit none

   integer, intent(in) :: i
   real, intent(in) :: biomass(Nmax,v_count,2), up(v_count,2)
   integer, intent(out) :: count_up
   integer :: count
   real :: M=0 !Tot mass in voxel

   M = sum(biomass(:,i,1) )
   call mass2cell_count(M, count)
   count_up = nint(count * up(i,1))
end

subroutine get_count_down(i,biomass,up,count_down)
   ! TODO Worth it?
   use input
   implicit none

   integer, intent(in) :: i
   real, intent(in) :: biomass(Nmax,v_count,2), up(v_count,2)
   integer, intent(out) :: count_down
   integer :: count
   real :: M=0 !Tot mass in voxel

   M = sum(biomass(:,i,1) )
   call mass2cell_count(M, count)
   count_down = nint(count * (1-up(i,1)))
end

subroutine get_count_particles(biomass_particle,eps_count_particle, count)
   ! TODO Worth it?
   use input
   implicit none

   integer, intent(in)  :: eps_count_particle
   real,    intent(in)  :: biomass_particle(Nmax)
   integer, intent(out) :: count
   real :: M=0
   M = sum(biomass_particle(:))

   call mass2cell_count(M,count)
   count = count + eps_count_particle
end

subroutine update_stochastics(i, c_q, biomass, up)
   use input !v_count, dt
   implicit none

   integer, intent(in)  :: i
   real,    intent(in)  :: c_q(v_count,2), biomass(Nmax,v_count,2)
   real,    intent(out) :: up(v_count,2)
   real :: d2u, u2d, rand
   integer :: count_up, count_down, count_d2u, count_u2d,j

   call get_count_down(i,biomass,up,count_down)
   call get_count_up(i,biomass,up,count_up)
   call probability_down2up(i,c_q,d2u)
   call probability_up2down(i,c_q,u2d)

   ! Down -> Up
   count_d2u = 0
   do j=1,count_down
      call random_number(rand)
      if (rand < dt * d2u) then
         count_d2u = count_d2u + 1
      end if
   end do
   ! Up -> Down
   count_u2d = 0
   do j=1,count_up
      call random_number(rand)
      if (rand < dt * u2d) then
         count_u2d = count_u2d + 1
      end if
   end do

   if (count_up > 0 .AND. up(i,1) > 0 ) then
      up(i,2) = up(i,1) + (count_d2u-count_u2d) * up(i,1)/real(count_up)
   else if (count_down > 0) then
      up(i,2) = (count_d2u - count_u2d) / real(count_down)
   end if
end

subroutine probability_down2up(i, c_q, prob)
   ! Calculates probability down to up
   use input !v_count, alpha, gamma
   implicit none
   integer, intent(in) :: i
   real, intent(in) :: c_q(v_count,2)
   real, intent(out) :: prob

   prob = alpha*c_q(i,1) / (1 + gamma *c_q(i,1))
end

subroutine probability_up2down(i, c_q, prob)
   ! Calculates probability down to up
   use input !v_count, beta, gamma
   implicit none
   integer, intent(in) :: i
   real, intent(in) :: c_q(v_count,2)
   real, intent(out) :: prob

   prob = beta / (1 + gamma *c_q(i,1))
end

subroutine update_production_s(i, c_s, biomass, prod_s)
   ! Calculate the new production of substrate and QSM
   ! in voxel i
   use input !Vmax, Ks, Nmax
   implicit none

   integer, intent(in) :: i
   real, intent(in) :: c_s(v_count,2), biomass(Nmax, v_count, 2)
   real, intent(out) :: prod_s(v_count)
   real :: M=0 ! Total mass
    
   M = sum(biomass(:,i,1))
   prod_s(i) = - Vmax* M * c_s(i,1) / (Ks + c_s(i,1))
end

subroutine update_concentration(i, prod, diff, c)
   ! Updates conentration (Forward Euler step)
   ! Input index, diffusion constant, voxel width, c before/after
   ! Output new concentration
   use input
   implicit none
   integer, intent(in) :: i
   real, dimension(v_count,2), intent(out) :: c
   real, intent(in) :: prod(v_count), diff
   integer :: j_list_neigh(6), j
   
   call get_index_neighbours(i, j_list_neigh)
   do j = 1,6 !For every neighbour
      if (.NOT. (j_list_neigh(j)  < 1 ) ) then
         c(i,2) = c(i,2) + dt * (c(j_list_neigh(j),1) - c(i,1)) *  diff/(v_width*v_width) + dt*prod(i)/(v_width**3)
      end if
   end do
end

subroutine append_biomass(mass, i, biomass)
   ! TODO check total particle count
   use input !Nmax
   implicit none

   integer, save, dimension(v_count) :: current=1 !Current index for voxel, to not overwrite data
   integer, intent(in) :: i
   real, intent(in) :: mass
   real, intent(out) :: biomass(Nmax,v_count,2)

   if (current(i) <= Nmax) then
      biomass(current(i), i, 2) = mass
      current(i) = current(i) + 1
   else
      print*, "Overfilled voxel biomass", i, current(i)
   end if
end

subroutine index2xyz(i, pos)
   ! Input index i
   ! Returns x y z in pos list
   use input
   implicit none
   integer, intent(in) :: i
   integer :: pos(3)
   integer :: x, y, z

   x = mod(i-1, v_size(1))
   y = mod( (i-1)/v_size(1), v_size(2))
   z = (i-1)/(v_size(1) * v_size(2))

   pos = (/x, y, z/)
end

subroutine xyz2index(pos, i)
   ! Input xyz in pos list
   ! Returns integer
   use input
   implicit none
   integer, intent(in) :: pos(3)
   integer, intent(out) :: i
   integer :: x,y,z
   x = pos(1)
   y = pos(2)
   z = pos(3)

   i = x + y * v_size(1) + z * v_size(1) * v_size(2) + 1
   
   if ((i < 1) .OR. (i > v_size(1)*v_size(2)*v_size(3))) then
      print*, "Error, i<=0 or i > v_count",i,x,y,z
      i = 0
   end if
end

subroutine get_index_neighbours(i, i_list_neigh)
   ! Input index
   ! Return list of indexes to neighbours
   use input
   implicit none
   integer, intent(in):: i
   integer, intent(out) :: i_list_neigh(6)
   integer :: pos(3), npos(3)
   i_list_neigh = (/0,0,0,0,0,0/)

   call index2xyz(i,pos) 

   ! x +- 1 (Boundary condition)
   npos = pos + (/1,0,0/)
   call continuous_boundary_condition(npos(1), v_size(1))
   call xyz2index(npos, i_list_neigh(1))

   npos = pos - (/1,0,0/)
   call continuous_boundary_condition(npos(1), v_size(1))
   call xyz2index(npos, i_list_neigh(2))

   ! y +- 1
   npos = pos + (/0,1,0/)
   call continuous_boundary_condition(npos(2), v_size(2))
   call xyz2index(npos, i_list_neigh(3))

   npos = pos - (/0,1,0/)
   call continuous_boundary_condition(npos(2), v_size(2))
   call xyz2index(npos, i_list_neigh(4))

   ! z +- 1
   npos = pos + (/0,0,1/)
   if (.NOT. (npos(3) >= v_size(3))) then
      call xyz2index(npos, i_list_neigh(5))
   end if

   npos = pos - (/0,0,1/)
   if (.NOT. (npos(3) < 0)) then
      call xyz2index(npos, i_list_neigh(6))
   end if
end

subroutine continuous_boundary_condition(pos_q, v_size_q)
   ! Input position in one direction with voxel size,
   ! Returns new position according to CBC
   implicit none
   integer, intent(out)::pos_q
   integer, intent(in)::v_size_q

   if (pos_q < 0) then
      pos_q = v_size_q-1
   else if (pos_q >= v_size_q) then
      pos_q = 0
   end if
end

subroutine mass2cell_count(mass, count)
   use input
   implicit none

   real, intent(in):: mass
   integer, intent(out) :: count

   count = ceiling(mass/avg_mass_cell)
end

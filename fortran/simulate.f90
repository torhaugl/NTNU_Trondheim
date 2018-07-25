module input
   implicit none

   ! Input file / Parameters
   real, parameter :: dt = 0.07/60.0 !time step (min)
   real, parameter :: final_time = 0.1*60.0 !minutes
   integer, parameter :: NumTrials = final_time / dt
   real, parameter :: c_bulk = 0.2 ! Concentration bulk substrate
   real, parameter :: v_width = 17.0 !Voxel length
   real, parameter :: diff_s = 40680.0, diff_q = 33300.0 !Diffusion subst./QSM
   integer, parameter :: v_size(3) = (/10,10,100/)
   integer, parameter :: v_count = v_size(1) * v_size(2) * v_size(3)

   real, parameter :: Vmax = 46e-3 !Maximum substrate uptake
   real, parameter :: Ks = 2.34e-3 ! Half-saturaton const (substrate uptake)

   real, parameter :: Zqd = 8.3, Zqu = 1230.0 !QSM production TODO Set values
   real, parameter :: Kq = 1 ! TODO Set value

   real, parameter :: Ymax = 0.444 ! TODO Set value
   real, parameter :: maintenance = 0.6e-3 !TODO Set value

   real, parameter :: avg_mass_cell = 420.0

   real, parameter :: density_cell = 290.0
   real, parameter :: max_vol_frac = 0.52
   real, parameter :: Mmax = 14700.0
   integer, parameter :: Nmax = floor(density_cell * max_vol_frac * v_width**3 / Mmax) ! Standard input: 50

end


program simulate

   ! Simulates time-steps
   use input
   implicit none

   ! Functions
   real, external :: avg

   ! Variables
   ! Concentration substrate (before and after timestep)
   real, dimension(v_count, 2) :: c_s, c_q, up
   real, dimension(Nmax, v_count, 2) :: biomass, epsmass

   real, dimension(v_count) :: prod_s, prod_q
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
   epsmass(:,:,:) = 0.0

   do i = 1, v_size(1)*v_size(2)
      call append_biomass(600.0,i, biomass)
      biomass(:,i,1) = biomass(:,i,2)
   end do

   !!! CODE
   print*,"CODE START"
   print*,"Number of steps:", NumTrials
   print*,

   call cpu_time(start)

   do n = 1,NumTrials
      do i = 1, v_count
         call update_mass(i, c_s, biomass)
         call update_production_s(i, c_s, biomass, prod_s) !Substrate
         call update_production_q(i, c_q, biomass, up, prod_q)
         call update_concentration(i, prod_s, diff_s, c_s) !Substrate
         call update_concentration(i, prod_q, diff_q, c_q) !QSM
      end do
      c_s(:,1) = c_s(:,2) ! Insert the newly calculated concentrations
      c_q(:,1) = c_q(:,2)
      biomass(:,:,1) = biomass(:,:,2)
      epsmass(:,:,1) = epsmass(:,:,2)
   end do

   call cpu_time(finish)

   print*,"CPU time(s):", finish-start
   print*,"Model time(min):", final_time
   print*,
   print*,"Top"   ,c_s(10000,1)  ,c_q(10000,1)
   print*,"Mid"   ,c_s(5000,1)   ,c_q(5000,1)
   print*,"Bot"   ,c_s(1,1)      ,c_q(1,1)
   print*,"AVG"   ,avg(c_s(:,1)) ,avg(c_q(:,1))
   print*,"Bio"   ,biomass(:4,1,1)

end program simulate

real function avg(arr)
   ! Must be length v_count
   use input
   implicit none
   real, intent(in) :: arr(v_count)
   avg = sum(arr) / size(arr)
   return
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
   real :: M=0 !Tot mass
   integer :: count_up, count_down

   M = sum(biomass(:,i,1) )
   call mass2cell_count(M*up(i,1), count_up)
   call mass2cell_count(M*(1-up(i,1)), count_down)

   prod_q(i) = Zqd*count_down + Zqu*count_up * c_q(i,1)/(Kq + c_q(i,1) )
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
   ! Input index, diffusion constant, vortex width, c before/after
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
   use input !Nmax
   implicit none

   integer, save, dimension(v_count) :: current=1 !Current index for voxel, to not overwrite data
   integer, intent(in) :: i
   real, intent(in) :: mass
   real, intent(out) :: biomass(Nmax,v_count,2)

   biomass(current(i), i, 2) = mass
   current(i) = current(i) + 1
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
   ! Input position in one direction with vortex size,
   ! Returns new position according to CBC
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
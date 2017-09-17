  ! A 1d profile model of river channel evolution

  Program 1dprofile
    use implicit none

    real :: transport, c, D, time, factor, max, totsed, delta_h, sum;
    integer :: bedrock, xc, y, lattice_size, check

    real, dimension (:), allocatable :: h, h_old
    
    !open(file="space1.txt")
    !open(file="time1.txt")

    c = 0.1
    D = 1.0
    lattice_size = 128 ! delta (D) is 1km, so this is a 128km basin
    bedrock = 16       ! bedrock-alluvial transition distance from divide

    allocate (h(lattice_size))
    allocate (h_old(lattice_size)

    xc = bedrock

    do ci = i, lattice_size
       if (i .gt. bedrock) then
          h(i) = 0.0
          h_old = 0.0

          ! plot the initial condition
          ! write(...)
       else
          h(i) = 1.0
          h_old(i) = 1.0
          ! write the value write(...)
       endif
    enddo

    factor = 1.0
    time = 0.0
    check = 1000
    
    do while (time .lt. 100000.0)
       if (time .gt. check) then
          sum = 0
          check += 10000
          do i=1, lattice_size
             sum += h(i)
             ! write result
          enddo

          ! write time sum/lattice-size
          ! write xc

          max = 0
          totsed = 0
          delta_h = c * factor * (1/float(lattice_size)) * (h(1) - h(2))
          if (abs(delta_h) .gt. max) then
             max = abs(delta_h)
          endif
          h(1) -= delta_h
          totsed = delta_h

          do i=2, lattice_size
             delta_h = c * factor * (i/float(lattice_size)) * h_old(i) - h_old(i+1)
             totsed += delta_h
             transport = D * factor * (i/float(lattice_size)) * h_old(i) - hold(i+1)
             if ((transport .gt. (totsed + 0.01)) .and. (i .le. bedrock)) then
                xc = i
                if (abs(delta_h) .gt. max) then
                   max = abs(delta_h)
                   h(i) -= delta_h
                else
                   delta_h = D * factor * ((i-1)/float(lattice_size)) * (h_old(i-1) - h_old(i)) - D * factor * (i/float(lattice_size) * (h_old(i) - h_old(i+1))
                   if (abs(delta_h) .gt. max) then
                      max = abs(delta_h)
                      h(i) -= delta_h
                   endif
                endif
             endif
          enddo

          if (max .lt. 0.001) then
             do i=1, lattice_size
                h_old(i) = h(i)
                time += factor
             enddo
          else
             h(i) = h_old(i)
             factor /= 3
          endif
          if (max .lt. 0.0001) then
             factor *= 3
          endif
          
                
       

module monte_carlo
    implicit none 
contains

    subroutine initial_setup(l, num_part, step)
        integer, intent(in) :: l, num_part
        integer, intent(out) :: step(l)
        integer :: i, k
        double precision :: rnd

        step = 0
        
        do i = 1, num_part
            call random_number(rnd)
            k = floor(rnd*dble(l)) + 1 
            step(k) = step(k) + 1
        end do

    end subroutine initial_setup

    subroutine random_move(l, step_in, step_out)
        integer, intent(in) :: l, step_in(l)
        integer, intent(out) :: step_out(l)
        double precision :: u, v
        integer :: k

        step_out = step_in

        call random_number(v)
        k = floor(dble(l)*v) + 1
        
        if (step_out(k)>0) then
            step_out(k) = step_out(k) - 1
            call random_number(u)
            if (u < 0.5) then
                if (k == l) then
                    step_out(1) = step_out(1) + 1
                else 
                    step_out(k+1) = step_out(k+1) + 1
                end if
            else 
                if (k == 1) then
                    step_out(l) = step_out(l) + 1
                else 
                    step_out(k-1) = step_out(k-1) + 1
                end if 
        
            end if 
        end if 

    end subroutine random_move

    subroutine random_move_ratio(l, step_in, step_out)
        integer, intent(in) :: l, step_in(l)
        integer, intent(out) :: step_out(l)
        double precision :: u, v
        integer :: k, i

        step_out = step_in
        do i = 0,l
            call random_number(v)
            k = floor(dble(l)*v) + 1
            
            if (step_out(k)>0) then
                step_out(k) = step_out(k) - 1
                call random_number(u)
                if (u < 0.5) then
                    if (k == l) then
                        step_out(1) = step_out(1) + 1
                    else 
                        step_out(k+1) = step_out(k+1) + 1
                    end if
                else 
                    if (k == 1) then
                        step_out(l) = step_out(l) + 1
                    else 
                        step_out(k-1) = step_out(k-1) + 1
                    end if 
            
                end if 
            exit
            end if 
        end do 

    end subroutine random_move_ratio
        

    subroutine random_number_change(l, step_in, step_out, adding)
        integer, intent(in) :: l, step_in(l)
        integer, intent(out) :: step_out(l)
        double precision :: u, v
        logical, intent(out) :: adding
        integer :: j, k, i, n_valid
        integer, allocatable :: valid_sites(:)

        n_valid = 0
        
        do i = 1, l
            if (step_in(i) > 0) then
                n_valid = n_valid + 1
            end if
        end do
        
        allocate(valid_sites(n_valid))
        valid_sites = 0
        n_valid = 0
        do i = 1, l
            if (step_in(i) > 0) then
                n_valid = n_valid + 1
                valid_sites(n_valid) = i
            end if
        end do
        step_out = step_in

        call random_number(v)
        if (v < 0.5) then
            adding = .TRUE. 
            call random_number(u)
            j = floor(dble(l)*u) + 1
            step_out(j) = step_out(j) + 1
        else 
            adding = .FALSE. 
            call random_number(u)
            k = valid_sites(floor(dble(n_valid)*u) + 1)
            step_out(k) = step_out(k) - 1
        end if

        deallocate(valid_sites)

    end subroutine random_number_change


    subroutine compute_energy(l, step_in, e)
        implicit none
        integer, intent(in) :: l
        integer, intent(in) :: step_in(l)
        double precision, intent(out) :: e
        integer :: j, n, k, r, nv, nh
        double precision :: hgs

        n = sum(step_in)
        k = n / l
        r = mod(n, l)

        nv = n !energia legami verticali
        nh = 0

        do j = 1, l-1
            nh = nh + min(step_in(j), step_in(j+1)) !energia legami orizzontali
        end do
        
        nh = nh + min(step_in(l), step_in(1)) !PBC
        e = -dble(nv + nh)
        hgs = -dble(n + l*k + max(0, r-1))
        e = e - hgs

    end subroutine compute_energy

    subroutine observables (l, step, kink_rho, step_energy, mean_len)
        integer, intent(in) :: l, step(l)
        double precision, intent(out) ::  kink_rho, step_energy, mean_len
        integer :: j
        double precision :: delta
        
        kink_rho = 0.0d0
        step_energy = 0.0d0
        mean_len = 0.0d0

        do j = 1, l
            mean_len = mean_len + (dble(step(j))/dble(l))
            if ( j /= l ) then
                delta = abs(step(j) - step(j+1))
                if (delta /= 0) then
                    kink_rho = kink_rho + (1.0d0/dble(l))
                end if 
                step_energy = step_energy + (delta/dble(l))
            else 
                delta = abs(step(j) - step(1))
                if (delta /= 0) then
                    kink_rho = kink_rho + (1.0d0/dble(l))
                end if 
                step_energy = step_energy + (delta/dble(l))
            end if 
        end do

    end subroutine observables 


    subroutine metropolis(l, T, step_old, step_new, energy, acc_move)
        integer, intent(in) :: l, step_old(l) 
        integer, intent(out) :: step_new(l)
        double precision, intent(in) :: T
        double precision, intent(out) :: acc_move, energy
        double precision :: e_trial, e_old, delta_e, w, rnd
        integer ::  step_trial(l)
    
        acc_move = 0
        call compute_energy(l, step_old, e_old)
        !print*, "CALCOLATA ENERGIA STEP VECCHIO"
        call random_move(l, step_old, step_trial)
        !print*,"STEP TRIAL", step_trial
        call compute_energy(l, step_trial, e_trial) 
        !print*, "CALCOLATA ENERGIA STEP TRIAL"
        delta_e = e_trial - e_old 
        
        w = exp((-1.0d0*delta_e)/T)
        !print*, "PROBABILITA", w
        if ( w >= 1.0d0 ) then
            step_new = step_trial
            acc_move = acc_move + 1
        else 
            call random_number(rnd)
            !print*, "NUMERO CASUALE", rnd
            if ( rnd <= w ) then 
                step_new = step_trial
                acc_move = acc_move + 1
            else 
                step_new = step_old
            end if
        end if 
        
        call compute_energy(l, step_new, energy)
        !print*, "CALCOLATA ENERGIA STEP AGGIORNATO"
    end subroutine metropolis

    subroutine metropolis_ratio(l, T, step_old, step_new, energy, acc_move)
        integer, intent(in) :: l, step_old(l) 
        integer, intent(out) :: step_new(l)
        double precision, intent(in) :: T
        double precision, intent(out) :: acc_move, energy
        double precision :: e_trial, e_old, delta_e, w, rnd, ratio
        integer ::  step_trial(l), m_old, m_new
    
        acc_move = 0
        call compute_energy(l, step_old, e_old)
        !print*, "CALCOLATA ENERGIA STEP VECCHIO"
        call random_move_ratio(l, step_old, step_trial)
        !print*,"STEP TRIAL", step_trial
        call compute_energy(l, step_trial, e_trial) 
        !print*, "CALCOLATA ENERGIA STEP TRIAL"
        delta_e = e_trial - e_old 
        
        m_old = count(step_old > 0)
        m_new = count(step_trial > 0)
        ratio = dble(m_old) / dble(m_new)

        w = exp((-1.0d0*delta_e)/T)*ratio
        !print*, "PROBABILITA", w
        if ( w >= 1.0d0 ) then
            step_new = step_trial
            acc_move = acc_move + 1
        else 
            call random_number(rnd)
            !print*, "NUMERO CASUALE", rnd
            if ( rnd <= w ) then 
                step_new = step_trial
                acc_move = acc_move + 1
            else 
                step_new = step_old
            end if
        end if 
        
        call compute_energy(l, step_new, energy)
        !print*, "CALCOLATA ENERGIA STEP AGGIORNATO"
    end subroutine metropolis_ratio

    subroutine metropolis_2T(l, T, step_old, step_new, energy, acc_move)
        integer, intent(in) :: l, step_old(l) 
        integer, intent(out) :: step_new(l)
        double precision, intent(in) :: T
        double precision, intent(out) :: acc_move, energy
        double precision :: e_trial, e_old, delta_e, w, rnd
        integer ::  step_trial(l)
    
        acc_move = 0
        call compute_energy(l, step_old, e_old)
        !print*, "CALCOLATA ENERGIA STEP VECCHIO"
        call random_move(l, step_old, step_trial)
        !print*,"STEP TRIAL", step_trial
        call compute_energy(l, step_trial, e_trial) 
        !print*, "CALCOLATA ENERGIA STEP TRIAL"
        delta_e = e_trial - e_old 

        w = exp((-1.0d0*delta_e)/(2*T))
        !print*, "PROBABILITA", w
        if ( w >= 1.0d0 ) then
            step_new = step_trial
            acc_move = acc_move + 1
        else 
            call random_number(rnd)
            !print*, "NUMERO CASUALE", rnd
            if ( rnd <= w ) then 
                step_new = step_trial
                acc_move = acc_move + 1
            else 
                step_new = step_old
            end if
        end if 
        
        call compute_energy(l, step_new, energy)
        !print*, "CALCOLATA ENERGIA STEP AGGIORNATO"
    end subroutine metropolis_2T


    subroutine metropolis_adding(l, T, step_old, step_new, energy, acc_move)
        integer, intent(in) :: l, step_old(l) 
        integer, intent(out) :: step_new(l)
        double precision, intent(in) :: T
        double precision, intent(out) :: acc_move, energy
        logical :: adding
        double precision :: e_trial, e_old, delta_e, w, rnd
        integer ::  step_trial(l), n_valid
    
        acc_move = 0

        call random_number_change(l, step_old, step_trial, adding)

        if (adding .eqv. .TRUE.) then
            step_new = step_trial
            acc_move = acc_move + 1
        else if (adding .eqv. .FALSE.) then
            
            n_valid = count(step_old > 0)
            call compute_energy(l, step_old, e_old)
            call compute_energy(l, step_trial, e_trial) 
            delta_e = e_trial - e_old 

            w = exp((-1.0d0*delta_e)/T)*(dble(n_valid)/dble(l))
            !print*, "PROBABILITA", w
            if ( w >= 1.0d0 ) then
                step_new = step_trial
                acc_move = acc_move + 1
            else 
                call random_number(rnd)
                !print*, "NUMERO CASUALE", rnd
                if ( rnd <= w ) then 
                    step_new = step_trial
                    acc_move = acc_move + 1
                else 
                    step_new = step_old
                end if
            end if 
        end if 
        call compute_energy(l, step_new, energy)
    end subroutine metropolis_adding

    subroutine n_fold(l, T, step_in, step_new, dt)
        integer, intent(in)  :: l, step_in(l)
        double precision, intent(in)  :: T
        integer, intent(out) :: step_new(l)
        double precision, intent(out) :: dt
        integer :: i, ip, im, ipp, imm, cls, chosen_class, chosen_idx
        integer :: delta_e_int
        double precision :: tot_prob, cum_prob, rnd, rnd2
        double precision :: class_rates(5)
        integer :: class_count(5)
        integer, allocatable :: class_from(:,:), class_to(:,:)

        class_rates(1) = 1.0d0
        class_rates(2) = 1.0d0
        class_rates(3) = 1.0d0
        class_rates(4) = exp(-1.0d0 / T)
        class_rates(5) = exp(-2.0d0 / T)

        allocate(class_from(5, 2*l))
        allocate(class_to  (5, 2*l))
        class_count = 0

        do i = 1, l
            if (step_in(i) < 1) cycle

            ip  = mod(i, l) + 1
            im  = mod(i-2+l, l) + 1
            ipp = mod(ip, l) + 1
            imm = mod(im-2+l, l) + 1
            
            delta_e_int = -((min(step_in(im),   step_in(i)-1) - min(step_in(im),  step_in(i)))  &
                        + (min(step_in(i)-1,  step_in(ip)+1) - min(step_in(i),  step_in(ip))) &
                        + (min(step_in(ip)+1, step_in(ipp))  - min(step_in(ip), step_in(ipp))))

            cls = delta_e_int + 3
            class_count(cls) = class_count(cls) + 1
            class_from(cls, class_count(cls)) = i
            class_to  (cls, class_count(cls)) = ip

            delta_e_int = -( (min(step_in(imm), step_in(im)+1) - min(step_in(imm), step_in(im))) &
                        + (min(step_in(im)+1, step_in(i)-1) - min(step_in(im), step_in(i)))  &
                        + (min(step_in(i)-1, step_in(ip)) - min(step_in(i), step_in(ip))))

            cls = delta_e_int + 3
            class_count(cls) = class_count(cls) + 1
            class_from(cls, class_count(cls)) = i
            class_to (cls, class_count(cls)) = im
        end do

        tot_prob = sum(class_count * class_rates)

        call random_number(rnd2)
        dt = -log(rnd2) / tot_prob

        call random_number(rnd)
        cum_prob = 0.0d0
        chosen_class = 5

        do cls = 1, 5
            cum_prob = cum_prob + class_count(cls) * class_rates(cls)
            if (rnd * tot_prob < cum_prob) then
                chosen_class = cls
                exit
            end if
        end do

        call random_number(rnd)
        chosen_idx = floor(rnd * dble(class_count(chosen_class))) + 1
        chosen_idx = min(chosen_idx, class_count(chosen_class))

        step_new = step_in
        step_new(class_from(chosen_class, chosen_idx)) = &
            step_new(class_from(chosen_class, chosen_idx)) - 1
        step_new(class_to  (chosen_class, chosen_idx)) = &
            step_new(class_to  (chosen_class, chosen_idx)) + 1

        deallocate(class_from, class_to)
    end subroutine n_fold

    subroutine autocorrelation_function(n_steps, taumax, kink_dens, step_en, energy, &
                                     aut_energy, aut_kink_den, aut_step_en)
        integer, intent(in) :: n_steps, taumax
        double precision, intent(in) :: kink_dens(n_steps), step_en(n_steps), energy(n_steps)
        double precision, intent(out) :: aut_energy(taumax), aut_kink_den(taumax), aut_step_en(taumax)
        integer :: i, tau, m
        double precision :: si_e, si2_e, sij_e
        double precision :: si_k, si2_k, sij_k
        double precision :: si_s, si2_s, sij_s

        aut_energy   = 0.0d0
        aut_kink_den = 0.0d0
        aut_step_en  = 0.0d0

        do tau = 0, taumax - 1
            m = n_steps - tau ! lunghezza del sotto-campione

            si_e  = 0.0d0 ; si2_e = 0.0d0 ; sij_e = 0.0d0
            si_k  = 0.0d0 ; si2_k = 0.0d0 ; sij_k = 0.0d0
            si_s  = 0.0d0 ; si2_s = 0.0d0 ; sij_s = 0.0d0

            do i = 1, m
                si_e  = si_e  + energy(i)
                si2_e = si2_e + energy(i)**2
                sij_e = sij_e + energy(i) * energy(i + tau)

                si_k  = si_k  + kink_dens(i)
                si2_k = si2_k + kink_dens(i)**2
                sij_k = sij_k + kink_dens(i) * kink_dens(i + tau)

                si_s  = si_s  + step_en(i)
                si2_s = si2_s + step_en(i)**2
                sij_s = sij_s + step_en(i) * step_en(i + tau)
            end do

            si_e  = si_e / dble(m) ; si2_e = si2_e / dble(m) ; sij_e = sij_e / dble(m)
            si_k  = si_k / dble(m) ; si2_k = si2_k / dble(m) ; sij_k = sij_k / dble(m)
            si_s  = si_s / dble(m) ; si2_s = si2_s / dble(m) ; sij_s = sij_s / dble(m)

            if (si2_e - si_e**2 > 1.0d-14) &
                aut_energy(tau+1) = (sij_e - si_e**2) / (si2_e - si_e**2)
            if (si2_k - si_k**2 > 1.0d-14) &
                aut_kink_den(tau+1) = (sij_k - si_k**2) / (si2_k - si_k**2)
            if (si2_s - si_s**2 > 1.0d-14) &
                aut_step_en(tau+1) = (sij_s - si_s**2) / (si2_s - si_s**2)
        end do

    end subroutine autocorrelation_function

    subroutine block_average(nsteps, blockstep, observable, err)

        integer, intent(in) :: nsteps, blockstep
        double precision, intent(in) :: observable(nsteps)
        double precision, intent(out) :: err
        integer :: nblock, i,j
        double precision :: mean_mean, variance
        double precision, allocatable :: mean(:)

        nblock = nsteps/blockstep
        if (nblock < 2) then
            err = 0.d0
            return
        end if

        allocate(mean(nblock))
        mean = 0.d0

        do i= 1, nblock
            do j= (i-1) * blockstep + 1, i*blockstep
                mean(i)= mean(i) + observable(j)
            end do
            mean(i)= mean(i)/dble(blockstep)
        end do

        mean_mean = sum(mean)/dble(nblock)
        variance = 0.d0

        do i=1,nblock
            variance = variance + (mean(i)-mean_mean)**2
        end do

        variance = variance/dble(nblock-1)
        err = sqrt(variance/dble(nblock))

        deallocate(mean)

    end subroutine block_average
end module monte_carlo

module main
    use monte_carlo
    implicit none
contains
    subroutine simulation(l, T, nsteps, nsteps_eq, num_part, move_type, seed_init, seed_mc, &
                          step_final, kink_dens, step_en, energy, mean_len, acceptance_prob, &
                          acceptance_prob1, aut_energy, aut_kink_den, aut_step_en, dt, nblk_pts, err_e, &
                          err_kink, err_step_e, blockstep_arr)
        integer, intent(in) :: l, nsteps, seed_init, seed_mc, num_part, nsteps_eq, nblk_pts
        integer, intent(out) :: step_final(l), blockstep_arr(nblk_pts)
        character(len = 32), intent(in) :: move_type 
        double precision, intent(in) :: T
        double precision, intent(out) :: kink_dens(nsteps), step_en(nsteps), energy(nsteps), dt(nsteps)
        double precision, intent(out) :: mean_len(nsteps), acceptance_prob, acceptance_prob1, &
                                        aut_energy(nsteps/80), aut_kink_den(nsteps/80), aut_step_en(nsteps/80)   
        double precision, intent(out) :: err_e(nblk_pts), err_kink(nblk_pts), err_step_e(nblk_pts)
        integer :: step_initial(l), i, j, k, idx, bmax, b
        double precision :: acc_move, sum_acc_move, sum_acc_move1, e_eq, dt_eq
        double precision :: total_time, mean_energy, mean_kink, mean_step_en
        integer :: seed_array(8)
        
        acceptance_prob = 0.d0
        acceptance_prob1 = 0.d0

        do i = 1, 8
            seed_array(i) = (seed_init + 3*i) * i * 2 - i + 273
        end do 
        call random_seed(put=seed_array)
        call initial_setup(l, num_part, step_initial)
        !print*, "step iniziale",step_initial
        
        do i = 1, 8
            seed_array(i) = (seed_mc + 10) * i * 2
        end do
        call random_seed(put=seed_array)

        if (move_type == "single_particle_move") then
            sum_acc_move = 0
            sum_acc_move1 = 0
            do j = 1, nsteps_eq
                call metropolis(l, T, step_initial, step_final, e_eq, acc_move)
            !    print*, "METROPOLIS FATTO"
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
            end do
            do i = 1,nsteps
            !    print*, "step iniziale check", step_initial
                call metropolis(l, T, step_initial, step_final, energy(i), acc_move)
                call observables(l, step_final, kink_dens(i), step_en(i), mean_len(i))
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
                sum_acc_move1 = sum_acc_move1 + acc_move
            !    print*, "step finale", step_final
            end do
            acceptance_prob = sum_acc_move/dble(nsteps + nsteps_eq)
            dt = 0
        else if (move_type == "single_particle_move_ratio") then
            sum_acc_move = 0
            sum_acc_move1 = 0
            do j = 1, nsteps_eq
                call metropolis_ratio(l, T, step_initial, step_final, e_eq, acc_move)
            !    print*, "METROPOLIS FATTO"
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
            end do
            do i = 1,nsteps
            !    print*, "step iniziale check", step_initial
                call metropolis_ratio(l, T, step_initial, step_final, energy(i), acc_move)
                call observables(l, step_final, kink_dens(i), step_en(i), mean_len(i))
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
                sum_acc_move1 = sum_acc_move1 + acc_move
            !    print*, "step finale", step_final
            end do
            acceptance_prob = sum_acc_move/dble(nsteps + nsteps_eq)
            dt = 0
        else if (move_type == "single_particle_move2T") then
            sum_acc_move = 0
            sum_acc_move1 = 0
            do j = 1, nsteps_eq
                call metropolis_2T(l, T, step_initial, step_final, e_eq, acc_move)
            !    print*, "METROPOLIS FATTO"
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
            end do
            do i = 1,nsteps
            !    print*, "step iniziale check", step_initial
                call metropolis_2T(l, T, step_initial, step_final, energy(i), acc_move)
                call observables(l, step_final, kink_dens(i), step_en(i), mean_len(i))
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
                sum_acc_move1 = sum_acc_move1 + acc_move
            !    print*, "step finale", step_final
            end do
            acceptance_prob = sum_acc_move/dble(nsteps + nsteps_eq)
            dt = 0
        else if (move_type == "nfold_move") then
            do j = 1,nsteps_eq
                call n_fold(l, T, step_initial, step_final, dt_eq)
            !   print*, "METROPOLIS FATTO"
                step_initial = step_final
            end do
            do i = 1,nsteps
            !   print*, step_initial
                call n_fold(l, T, step_initial, step_final, dt(i))
                call compute_energy(l, step_final, energy(i))
                call observables(l, step_final, kink_dens(i), step_en(i), mean_len(i))
                step_initial = step_final
            end do
            acceptance_prob = 1.0
            total_time = sum(dt)
            mean_energy = sum(energy * dt) / total_time
            mean_kink = sum(kink_dens * dt) / total_time
            mean_step_en = sum(step_en * dt) / total_time
        else if (move_type == "changing_number_particle") then
            sum_acc_move = 0
            sum_acc_move1 = 0
            do j = 1,nsteps_eq
                call metropolis_adding(l, T, step_initial, step_final, e_eq, acc_move)
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
            end do 
            do i = 1,nsteps
                call metropolis_adding(l, T, step_initial, step_final, energy(i), acc_move)
                call observables(l, step_final, kink_dens(i), step_en(i), mean_len(i))
                step_initial = step_final
                sum_acc_move = sum_acc_move + acc_move
                sum_acc_move1 = sum_acc_move1 + acc_move
            end do
            acceptance_prob = sum_acc_move/dble(nsteps + nsteps_eq)
            acceptance_prob1 = sum_acc_move1/dble(nsteps)
            dt = 0
        else 
            print*, "Selezionare tipo di mossa: single_particle_move / single_particle_move_ratio / "// &
            "single_particle_move2T / nfold_move / changing_number_particle"
            stop 1
        end if
        
        call autocorrelation_function(nsteps, nsteps/80, kink_dens, step_en, &
                                        energy, aut_energy, aut_kink_den, aut_step_en)
        
        bmax = max(1, nsteps/10) 
        do idx = 1, nblk_pts
            b = 1 + (idx-1)*(bmax-1)/max(1, nblk_pts-1)
            blockstep_arr(idx) = b
            call block_average(nsteps, b, energy, err_e(idx))
            call block_average(nsteps, b, kink_dens, err_kink(idx))
            call block_average(nsteps, b, step_en, err_step_e(idx))
        end do

    end subroutine simulation
end module main
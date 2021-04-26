        program calc_flip
        implicit real*8(a-h,o-z)
        parameter (npermut=720)
        parameter (noxy=6)
        parameter (nmax=100000)
        parameter (ndim=54)
        parameter (nwavetot=100)
        parameter (nmaxtot=nwavetot*nmax)
        parameter (natoms=18)
        parameter (nh=12)
        parameter (nbin=75)
        dimension nwater(noxy,npermut),iperm(nmaxtot),n(nwavetot),
     1  n0(nwavetot),time(nwavetot),psips(ndim,nmax),roox1_1(nmax),
     1  rooy1_1(nmax),rooz1_1(nmax),roox2_1(nmax),rooy2_1(nmax),
     1  rooz2_1(nmax),roox3_1(nmax),rooy3_1(nmax),rooz3_1(nmax),
     1  x(3,natoms),cross1(3,nmax),roox1_2(nmax),rooy1_2(nmax),
     1  rooz1_2(nmax),roox2_2(nmax),rooy2_2(nmax),rooz2_2(nmax),
     1  roox3_2(nmax),rooy3_2(nmax),rooz3_2(nmax),cross2(3,nmax),
     1  cross_norm1(nmax),cross_norm2(nmax),roh1(nmax),roh2(nmax),
     1  weight(nmax),r_coord1(3,nmax),r_coord2(3,nmax),nfh1(nmax),
     1  nfh2(nmax),rval1(nbin**2),box(nbin**2),dot1(nmax),
     1  coord(3,natoms,nmax),dot_coord(nmax),nfh_coord(nmax),
     1  roh_coord(nmax),tot_box(nbin**2),dot2(nmax),rval2(nbin**2),
     1  trash(ndim,nmax),weight_trash(nmax),box_max_sim(nwavetot),
     1  quad1(nbin**2),quad2(nbin**2),quad3(nbin**2),quad4(nbin**2),
     1  rval1tot(nbin**2),rval2tot(nbin**2)
        open(unit=9,file='wf-50k_walkers.dat',status='old',
     1  form='unformatted')
        open(unit=11,file='tot_weights_50k_walkers.dat',status='old')
        open(unit=12,file='flip_2d.dat',status='unknown')
C		This calculations takes the walkers and the descendant weights of each walker 
C		from a DMC calculation and generates projections of 4 outer OH bonds flipping
C		from an either an UB, DB,UF, DF configurations. 
C		Inputs: wfns
C		number of wave functions
C		number of walkers at time t, number of walkers at t=0, time t
C		Cartesian coordinate of all the atoms in the walkers from simulations 
C		in units of bohr (OHH,OHH,OHH,OHH,OHH,OHH)
C		weights 
C		number of walkers at time t
C		descendant weights of walkers
C		Outputs:
C		Projections of the flipping coordinates of the outer water molecule
        do i = 1,nbin**2
            tot_box(i) = 0.d0
        enddo
        read(9) nwave
        ip = 0
        it = 0
        do k = 1,nwave
            read(9) n(k),n0(k),time(k)
            do i = 1,n(k)
                ip = ip + 1
                read(9) (psips(j,i),j=1,ndim)
            enddo
            read(11,*) n1
            do i = 1,n1
                it = it + 1
                read(11,*) weight(i)
c                weight_tot(it) = weight(i)
            enddo
c            ntot = ip
c            tot_weight = 0.d0
c            do i = 1,ntot
c                tot_weight = tot_weight + weight_tot(i)
c            enddo
            do i = 1,n(k)
                is = 0
                do j = 1,natoms
                    do l = 1,3
                        is = is + 1
                        coord(l,j,i) = psips(is,i)*0.52917721067d0
                    enddo
                enddo
            enddo
c       calculate OOO plane for 1 of 2
            do i = 1,n(k)
                roox1_1(i) = psips(19,i)-psips(1,i)
                rooy1_1(i) = psips(20,i)-psips(2,i)
                rooz1_1(i) = psips(21,i)-psips(3,i)
                roox2_1(i) = psips(37,i)-psips(19,i)
                rooy2_1(i) = psips(38,i)-psips(20,i)
                rooz2_1(i) = psips(39,i)-psips(21,i)
                roox3_1(i) = psips(37,i)-psips(1,i)
                rooy3_1(i) = psips(38,i)-psips(2,i)
                rooz3_1(i) = psips(39,i)-psips(3,i)
                roox1_2(i) = psips(28,i)-psips(10,i)
                rooy1_2(i) = psips(29,i)-psips(11,i)
                rooz1_2(i) = psips(30,i)-psips(12,i)
                roox2_2(i) = psips(46,i)-psips(28,i)
                rooy2_2(i) = psips(47,i)-psips(29,i)
                rooz2_2(i) = psips(48,i)-psips(30,i)
                roox3_2(i) = psips(46,i)-psips(10,i)
                rooy3_2(i) = psips(47,i)-psips(11,i)
                rooz3_2(i) = psips(48,i)-psips(12,i)
            enddo
            do i = 1,n(k)
                cross_norm1(i) = 0.d0
                cross_norm2(i) = 0.d0
            enddo
            do i = 1,n(k)
                cross1(1,i) = rooy1_1(i)*rooz2_1(i)-(rooz1_1(i)*
     1          rooy2_1(i))
                cross1(2,i)=-(roox1_1(i)*rooz2_1(i)-(rooz1_1(i)*
     1          roox2_1(i)))
                cross1(3,i) = roox1_1(i)*rooy2_1(i)-(rooy1_1(i)*
     1          roox2_1(i))
                cross2(1,i) = rooy1_2(i)*rooz2_2(i)-(rooz1_2(i)*
     1          rooy2_2(i))
                cross2(2,i)=-(roox1_2(i)*rooz2_2(i)-(rooz1_2(i)*
     1          roox2_2(i)))
                cross2(3,i) = roox1_2(i)*rooy2_2(i)-(rooy1_2(i)*
     1          roox2_2(i))
            enddo
            do i = 1,n(k)
                do j = 1,3
                    cross_norm1(i) = cross_norm1(i)+(cross1(j,i)**2)
                    cross_norm2(i) = cross_norm2(i)+(cross2(j,i)**2)
                enddo
            enddo
            do i = 1,n(k)
                do j = 1,3
                    cross1(j,i) = cross1(j,i)/sqrt(cross_norm1(i))
                    cross2(j,i) = cross2(j,i)/sqrt(cross_norm2(i))
                enddo
            enddo
            call calc_roh_free(n(k),psips,5,roh1,r_coord1,nfh1)
            call calc_roh_free(n(k),psips,6,roh2,r_coord2,nfh2)
            call calc_dot_product(n(k),cross1,r_coord1,dot1)
            call calc_dot_product(n(k),cross2,r_coord2,dot2)
c            do i = 1,n(k)
c                if ((dot1(i)*0.52917721067d0.gt.0.d0).and.(dot2(i)*
c     1          0.52917721067d0.gt.0.d0)) then
c                    print *, natoms
c                    print *,dot1(i)*0.52917721067,dot2(i)*0.52917721067,
c     1              i,n(k),k
c                    do j = 1,natoms
c                        if (mod(j,3).eq.1) then
c                            print *, 'O',(coord(l,j,i),l=1,3)
c                        else
c                            print *, 'H',(coord(l,j,i),l=1,3)
c                        endif
c                    enddo
c                    print *, ' '
c                    read(*,*)
c                endif
c            enddo
            call descendant_weight_1(dot1(:),dot2(:),weight(:),n(k),
     1      rval1(:),rval2(:),box(:),dx1,dx2,nbin)
            do i =1,nbin**2
                tot_box(i) = tot_box(i) + box(i)
            enddo
            do i = 1,nbin**2
                box(i) = 0.d0
            enddo
        enddo
        tot_box = tot_box/dfloat(nwave)
        box_max = maxval(tot_box,1)
        do i = 1,nbin**2
            tot_box(i) = tot_box(i)/box_max
        enddo
        it = 0
        do i = 1,nbin
            do j = 1,nbin
                it = it + 1
                rval1tot(it) = rval1(i)
                rval2tot(it) = rval2(j)
                write(12,*) rval1(i),rval2(j),tot_box(it)
                tot_val = tot_val + tot_box(it)
            enddo
        enddo
        it1 = 0
        it2 = 0
        it3 = 0
        it4 = 0
        tot_quad1 = 0.d0
        tot_quad2 = 0.d0
        tot_quad3 = 0.d0
        tot_quad4 = 0.d0
        do i = 1,nbin**2
            if ((rval1tot(i).lt.0.d0).and.(rval2tot(i).lt.0.d0)) then
                it1 = it1 + 1
                quad1(it1) = tot_box(i)
                tot_quad1 = tot_quad1 + quad1(it1)
            else if((rval1tot(i).gt.0.d0).and.(rval2tot(i).lt.0.d0))then
                it2 = it2 + 1
                quad2(it2) = tot_box(i)
                tot_quad2 = tot_quad2 + quad2(it2)
            else if((rval1tot(i).lt.0.d0).and.(rval2tot(i).gt.0.d0))then
                it3 = it3 + 1
                quad3(it3) = tot_box(i)
                tot_quad3 = tot_quad3 + quad3(it3)
            else if((rval1tot(i).gt.0.d0).and.(rval2tot(i).gt.0.d0))then
                it4 = it4 + 1
                quad4(it4) = tot_box(i)
                tot_quad4 = tot_quad4 + quad4(it4)
            endif
        enddo
        avg_quad1 = tot_quad1/tot_val
        avg_quad2 = tot_quad2/tot_val
        avg_quad3 = tot_quad3/tot_val
        avg_quad4 = tot_quad4/tot_val
		print *, 'DB','UB','DF','UF'
        print *, avg_quad1,avg_quad2,avg_quad3,avg_quad4
        end program




        subroutine calc_roh_free(n,psips,key,r_oh,r_coord,nh)
        implicit real*8(a-h,o-z)
        parameter (nmax=100000)
        parameter (ndim=54)
        dimension psips(ndim,nmax),sec_roh1(nmax),sec_roh2(nmax),
     1  sec_roh3(nmax),sec_roh4(nmax),sec_roh_1(nmax),sec_roh_2(nmax),
     1  r_oh(nmax),r_coord(3,nmax),nh(nmax),test(4,nmax)
C		Calculation defines two planes and locates outer water molecule
C		with an free OH flipping back and forth.
C		In these calculations, assume that the 5th and 6th water molecule
C		contains the free OH
C		Inputs:
C		n = number of walkers
C		psips = 3*natoms, nwalker array of water clusters in units of bohr 
C		key = targeting specific water molecule (number 5-6)
C		Outputs:
C		r_oh = array of nwalker roh bonds that contain the free water molecule
C		r_coord = 3,nwalker Cartesian coordinates that contain the free
C		water molecule
C		nh = array identifying OH bond (either 1 or 2)
        do i = 1,n
            if (key.eq.5) then
                sec_roh1(i) = sqrt((psips(40,i)-psips(1,i))**2+
     1          (psips(41,i)-psips(2,i))**2+(psips(42,i)-psips(3,i))**2)
                sec_roh2(i) = sqrt((psips(40,i)-psips(19,i))**2+(
     1          psips(41,i)-psips(20,i))**2+(psips(42,i)-
     1          psips(21,i))**2)
                sec_roh3(i) = sqrt((psips(43,i)-psips(1,i))**2+(
     1          psips(44,i)-psips(2,i))**2+(psips(43,i)-psips(3,i))**2)
                sec_roh4(i) = sqrt((psips(43,i)-psips(19,i))**2+(
     1          psips(44,i)-psips(20,i))**2+(psips(45,i)-
     1          psips(21,i))**2)
            else if (key.eq.6) then
                sec_roh1(i) = sqrt((psips(49,i)-psips(10,i))**2+
     1          (psips(50,i)-psips(11,i))**2+(psips(51,i)-
     1          psips(12,i))**2)
                sec_roh2(i) = sqrt((psips(49,i)-psips(28,i))**2+(
     1          psips(50,i)-psips(29,i))**2+(psips(51,i)-
     1          psips(30,i))**2)
                sec_roh3(i) = sqrt((psips(52,i)-psips(10,i))**2+(
     1          psips(53,i)-psips(11,i))**2+(psips(54,i)-
     1          psips(12,i))**2)
                sec_roh4(i) = sqrt((psips(52,i)-psips(28,i))**2+(
     1          psips(53,i)-psips(29,i))**2+(psips(54,i)-
     1          psips(30,i))**2)
            endif
        enddo
        do i = 1,n
            if (sec_roh1(i).lt.sec_roh2(i)) then
                sec_roh_1(i) = sec_roh1(i)
            else
                sec_roh_1(i) = sec_roh2(i)
            endif
            if (sec_roh3(i).lt.sec_roh4(i)) then
                sec_roh_2(i) = sec_roh3(i)
            else
                sec_roh_2(i) = sec_roh4(i)
            endif
            if (key.eq.5) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(43,i)-psips(37,i))**2+(
     1              psips(44,i)-psips(38,i))**2+(psips(45,i)-
     1              psips(39,i))**2)
                    r_coord(1,i) = psips(43,i)-psips(37,i)
                    r_coord(2,i) = psips(44,i)-psips(38,i)
                    r_coord(3,i) = psips(45,i)-psips(39,i)
                    nh(i) = 2
                else
                    r_oh(i) = sqrt((psips(40,i)-psips(37,i))**2+(
     1              psips(41,i)-psips(38,i))**2+(psips(42,i)-
     1              psips(39,i))**2)
                    r_coord(1,i) = psips(40,i)-psips(37,i)
                    r_coord(2,i) = psips(41,i)-psips(38,i)
                    r_coord(3,i) = psips(42,i)-psips(39,i)
                    nh(i) = 1
                endif
            else if (key.eq.6) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(52,i)-psips(46,i))**2+(
     1              psips(53,i)-psips(47,i))**2+(psips(54,i)-
     1              psips(48,i))**2)
                    r_coord(1,i) = psips(52,i)-psips(46,i)
                    r_coord(2,i) = psips(53,i)-psips(47,i)
                    r_coord(3,i) = psips(54,i)-psips(48,i)
                    nh(i) = 2
                else
                    r_oh(i) = sqrt((psips(49,i)-psips(46,i))**2+(
     1              psips(50,i)-psips(47,i))**2+(psips(51,i)-
     1              psips(48,i))**2)
                    r_coord(1,i) = psips(49,i)-psips(46,i)
                    r_coord(2,i) = psips(50,i)-psips(47,i)
                    r_coord(3,i) = psips(51,i)-psips(48,i)
                    nh(i) = 1
                endif
            endif
        enddo
        return
        end subroutine

        subroutine calc_dot_product(n,vec1,vec2,dot)
        implicit real*8(a-h,o-z)
        parameter (nmax=100000)
        dimension vec1(3,nmax),vec2(3,nmax),dot(nmax)
C		Calculates the dot product of 2 vectors
C		Inputs:
C		n = number of vectors considered
C		vec1 = array of first vector
C		vec2 = array of second vector
C		Outputs:
C		dot = dot product of the two vectors
        do i = 1,n
            dot(i) = 0.d0
        enddo
        do i = 1,n
            do j = 1,3
                dot(i) = dot(i) + vec1(j,i)*vec2(j,i)
            enddo
        enddo
        return
        end subroutine

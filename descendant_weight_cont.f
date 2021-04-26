        subroutine descendant_weight_1(x,y,fam,n,xval,yval,
     1  standard_box,dx1,dx2,nbin)
        implicit real*8 (a-h,o-z)
C       this program will take my coordinates and generate a 2D
C		histogram based on the coordinates of interest.
C		Inputs: 
C		x = array of first coordinate of interest
C		y = array of second coordinate of interest
C		fam = array of weights calculated via descendant weighting
C		n = number of walkers
C		nbin = number of bins for histogram
C		Outputs:
C		xval = equally spaced coordinates along the x-axis
C		yval = equally spaced coordinates along the y-axis
C		dx1 = dx for x axis
C		dx2 = dx for y axis
C		standard_box = values of histogram based on xval and yval. 
C		Values are normalized such that total area under the curve of
C		\Psi^2 is 1.
        parameter (nmax=1000000)
        parameter (ndim=1)
        dimension x(nmax),y(nmax),nfam(nmax),xval(nbin),box(nbin,nbin),
     1  box_squared(nbin*nbin),standard_box(nbin*nbin),x_new(nmax),
     1  fam(nmax),psi_sq(nbin),box_stretch(nbin*nbin),
     1  yval(nbin),box_norm(nbin)
        pi = dacos(-1.d0)
        do i = 1,nbin
            do j = 1,nbin
                box(i,j) = 0.d0
            enddo
        enddo
        descend = 0.d0
        do i =1,n
            descend = descend + fam(i)
        enddo
        xmin = -3.d0
        xmax = 3.d0
        ymin = -3.d0
        ymax = 3.d0
        call hist_weight(x,y,nbin,n,box,xval,xmin,xmax,yval,ymin,ymax,
     1  fam,descend,dx1,dx2)
        ip = 0
        box_total = 0.d0
        do i = 1,nbin
            do j = 1,nbin
                ip = ip + 1
                box_stretch(ip) = box(i,j)
            enddo
        enddo
        do i = 1,nbin**2
            box_squared(i) = box_stretch(i)**2
            box_total = box_total + box_stretch(i)*dx1*dx2
        enddo
        do i = 1,nbin**2
            standard_box(i) = box_stretch(i)/box_total
        enddo
        tot_box = 0.d0
        do i = 1,nbin**2
            tot_box = tot_box + standard_box(i)*dx1*dx2
        enddo
c        print *, tot_box
        return
        end subroutine
            

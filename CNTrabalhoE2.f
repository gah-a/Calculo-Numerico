        implicit none
        double precision Int, a, b, fa, fb, e, C, pi, t(0:10), dt(1:10)
        double precision x, r(0:10)
        integer i, j
        
        open(10, file='dados2.txt', status='unknown')

        pi = acos(-1.d0)
        x = 1.d0
        e = 0.2d0
        t(0)=0.d0

            do i=1,10

                a = (i-1.d0)*pi/5.d0
                b = i*pi/5.d0

                fa = 1.d0/(1.d0+e*cos(a))**2.d0
                fb = 1.d0/(1.d0+e*cos(b))**2.d0

                Int = ((fb + fa)*(pi/5.d0))/2.d0

                C = 2.d0*pi*((1.d0-e**2.d0)**(-(3.d0/2.d0)))

                dt(i) = Int/C
                
                t(i) = dt(i) + t(i-1)
                
                r(i) = (x*(1.d0 - e**2.d0))/(1 + e*cos(a))

                write(*,*) 'O tempo', i, '=', t(i), 'valor de r', r(i)
                write(10,*) t(i), r(i)

            end do

        close(10)

        stop 
        end

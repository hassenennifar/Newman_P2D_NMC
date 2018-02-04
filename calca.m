function [a1, sumpa] = calca(k,ts, css, P)
% subroutine to calculate diffusion in solid particles
% of active material, based on a superposition integral.

a1 = zeros(P.nj,1);
sumpa = zeros(P.nj,1);

for j = 1:P.nj
    
    if j > P.bnd_sep_neg && j < P.bnd_pos_sep
        a1(j) = 0;
        sumpa(j) = 0;
    else

        if j <= P.bnd_sep_neg
            Rp = P.Rp_neg;
            Ds = P.Ds_neg;
        elseif j >= P.bnd_pos_sep
            Rp = P.Rp_pos;
            Ds = P.Ds_pos;
        end

        istart = 1;
        if ts(istart+1)-ts(istart) == 0
            istart = 2;
        end

        ai = zeros(k-1-istart,1);
        sumpa(j) = 0;

        s = pi^2/6;

        for i = istart:k-1
            tau(1) = Ds*(ts(k)-ts(i))/Rp^2;
            tau(2) = Ds*(ts(k)-ts(i+1))/Rp^2;

            an = zeros(2,1);
            for m = 1:2
                if tau(m) > 0.06
                    for n = 1:5
                        an(m) = an(m) + exp(-n^2*pi^2*tau(m))/n^2;
                    end
                    an(m) = 2/pi^2*(s-an(m));
                else
                    if tau(m) <= 0
                        an(m) = 0;
                    else
                        for n = 1:3
                            an(m) = an(m) + exp(-n^2/tau(m)) - n*sqrt(pi/tau(m))*erfc(n/sqrt(tau(m)));
                        end
                        an(m) = -tau(m)+2*sqrt(tau(m)/pi)*(1+2*an(m));
                    end
                end
            end

            an = an(1)-an(2);
            if i < k-1
                ai(k-1-i) = an;
                sumpa(j) = sumpa(j) + (css(j,i+1)-css(j,i))*ai(k-1-i)/(ts(i+1)-ts(i));
            else
                a1(j) = an;
            end
        end
    end
end

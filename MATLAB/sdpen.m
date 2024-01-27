%============================================================================================
%    SDPEN - A Sequential Penalty Derivative-free Method for
%    Nonlinear Constrained Optimization
%    Copyright (C) 2011  G.Liuzzi, S.Lucidi, M.Sciandrone
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    G. Liuzzi, S. Lucidi, M. Sciandrone. Sequential Penalty Derivative-free Methods for
%    Nonlinear Constrained Optimization, SIAM J. on Optimization, 20(5): 2614-2635 (2010)
%
%============================================================================================
function [xout, fout, ni, nf, eps, istop] = sdpen(foc,x,bl,bu,alfa_stop,nf_max,iprint,epsin)
    eps     = epsin;
    n       = length(bl);
    m       = length(eps);
    num_fal = 0;
    istop   = 0;
    nf      = 0;
    ni      = 0;
    fstop   = zeros(n+1,1);
    d       = ones(n,1);
    f       = funct(x);
    nf      = nf+1;
    alfa_d  = zeros(n,1);
    i_corr  = 1;
    z       = x;
    fz      = f;
    fstop(i_corr) = f;
    i_corr_fall   = 0;
    
    for i = 1:n
        alfa_d(i) = max(1.e-3,min(1.0,abs(x(i))));
        if (iprint >= 2)
            fprintf(' alfainiz(%d)=%e\n',i,alfa_d(i));
        end
    end
    alfa_max = max(alfa_d);
    
    if (iprint >= 2)
        fprintf(' ----------------------------------\n');
        fprintf(' finiz = %e\n',f);
        for i = 1:n
            fprintf(' xiniz(%d) = %e\n',i,x(i));
        end
    end

    while (true)

        if (iprint >= 0)
            fprintf(' ni = %d  nf = %d  f = %e  alfamax = %e\n',ni,nf,f,alfa_max);
        end
        if (iprint >= 2)
            for i = 1:n
                fprintf(' x(%d) = %e\n',i,x(i));
            end
        end
        
        %chiama la linesearch
        alfa = linesearch();
        
        if (abs(alfa) >= 1.e-12)
            x(i_corr)     = x(i_corr) + alfa*d(i_corr);
            f             = fz;
            fstop(i_corr) = f;
            num_fal       = 0;
        else
            if (i_corr_fall < 2)
                fstop(i_corr) = fz;
                num_fal       = num_fal+1;
            end
        end

        ni = ni+1;
        z(i_corr) = x(i_corr);

        if (i_corr < n)
            i_corr = i_corr+1;
        else
            i_corr = 1;
        end
        
        %call stop(n,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max)
        istop = stop();
        
        if (istop >= 1)
            break
        end

        %------------------------------------------------
        % Aggiornamento parametro di smoothing eps
        %------------------------------------------------
        cambio_eps = false;
        maxeps     = max(eps);
        for i = 1:m
            if (eps(i) == maxeps)
                if (eps(i) > 1.0e-2*sqrt(alfa_max))
                    if (iprint >= 1)
                        fprintf('**************************************\n');
                        fprintf('*********** aggiorno eps *************\n');
                    end
                    eps(i)     = min(1.e-2*eps(i), 1.0e-1*sqrt(alfa_max));
                    cambio_eps = true;
                    f          = funct(x);
                end
            end
        end
        if (cambio_eps)
            for i =1:n 
                alfa_d(i) = max(1.e-3,min(1.0,abs(x(i))));
            end
        end
        
    end
    
    if (iprint >= -1)
        fprintf(' ni = %d  nf = %d  f = %e  alfamax = %e\n',ni,nf,f,alfa_max);
    end
    fout = f;
    xout = x;
    return
    
    function fpen = funct(x)
        [fob,con] = feval(foc,x);
        
        fmax = 0.0;

        for i = 1:m
            fmax = fmax + (max(0.0,con(i))^1.1)/eps(i);
        end

        fpen = fob + fmax;
    end

    function istop = stop()

        istop    = 0;
        alfa_max = max(alfa_d);

        if (ni >= (n+1))
            ffm = f;
            for i = 1:n
                ffm = ffm+fstop(i);
            end
            ffm    = ffm/(n+1);
            ffstop = (f-ffm)*(f-ffm);
            for i = 1:n
                ffstop = ffstop+(fstop(i)-ffm)*(fstop(i)-ffm);
            end
            ffstop = sqrt(ffstop/(n+1));
        end

        if (alfa_max <= alfa_stop)
            istop = 1;
        end
        if (nf > nf_max)
            istop = 1;
        end

        return
    end

    function alfa = linesearch()

        gamma       = 1.e-6;
        delta       = 0.5;
        delta1      = 0.5;
        i_corr_fall = 0;
        ifront      = 0;
        j           = i_corr;

        if (iprint >= 1)
            fprintf(' j = %d    d(j) = %e\n',j,d(j));
        end

        if(abs(alfa_d(j)) <= 1.e-3*min(1.0,alfa_max))
            alfa = 0.0;
            if (iprint >= 1)
                fprintf('  alfa piccolo\n');
                fprintf(' alfa_d(j) = %e    alfamax = %e\n',alfa_d(j),alfa_max);
            end
            return
        end

        for ielle = 1:2

            if (d(j) > 0.0)
                if ((alfa_d(j)-(bu(j)-x(j))) < -1.e-6)
                    alfa   = max(1.e-24,alfa_d(j));
                else
                    alfa   = bu(j)-x(j);
                    ifront = 1;
                end
            else
                if ((alfa_d(j)-(x(j)-bl(j))) < -1.e-6)
                    alfa   = max(1.e-24,alfa_d(j));
                else
                    alfa   = x(j)-bl(j);
                    ifront = 1;
                end
            end

            if (abs(alfa) <= 1.e-3*min(1.0,alfa_max))
                d(j)        = -d(j);
                i_corr_fall = i_corr_fall+1;
                alfa        = 0.0;
                ifront      = 0;

                if (iprint >= 1)
                    fprintf(' direzione opposta per alfa piccolo\n');
                    fprintf(' j = %d    d(j) = %e\n',j,d(j));
                    fprintf(' alfa = %e    alfamax = %e\n',alfa,alfa_max);
                end

                continue
            end

            alfaex = alfa;
            z(j)   = x(j)+alfa*d(j);
            fz     = funct(z);
            nf     = nf+1;
            fpar   = f-gamma*alfa*alfa;


            if (iprint >= 1)
                fprintf(' fz = %e   alfa = %e\n',fz,alfa);
            end
            if (iprint >= 2)
                for i = 1:n
                    fprintf(' z(%d) = %e\n',i,z(i));
                end
            end

            if (fz < fpar)
                while (true)
                    if ((ifront == 1) || (num_fal > n-1))
                        alfa_d(j) = delta*alfa;
                        return
                    end

                    if (d(j) > 0.0)
                        if ((alfa/delta1-(bu(j)-x(j))) < -1.d-6)
                            alfaex = alfa/delta1;
                        else
                            alfaex = bu(j)-x(j);
                            ifront = 1;
                            if (iprint >= 1)
                                fprintf(' punto espan. sulla front.\n');
                            end
                        end
                    else
                        if ((alfa/delta1-(x(j)-bl(j))) < -1.d-6)
                            alfaex = alfa/delta1;
                        else
                            alfaex = x(j)-bl(j);
                            ifront = 1;
                            if (iprint >= 1)
                                fprintf(' punto espan. sulla front.\n');
                            end
                        end
                    end

                    z(j)    = x(j)+alfaex*d(j);
                    fzdelta = funct(z);
                    nf      = nf+1;
                    fpar    = f-gamma*alfaex*alfaex;

                    if (iprint >= 1)
                        fprintf(' fzex = %e  alfaex = %e\n',fzdelta,alfaex  );
                    end
                    if (iprint >= 2)
                        for i = 1:n
                            fprintf(' z(%d) = %e\n',i,z(i));
                        end
                    end

                    if (fzdelta < fpar)
                        fz        = fzdelta;
                        alfa      = alfaex;
                    else
                        alfa_d(j) = delta*alfa;
                        return
                    end
                end
            else

                d(j)   = -d(j);
                ifront = 0;

                if (iprint >= 1)
                    fprintf(' direzione opposta\n');
                    fprintf(' j = %d    d(j) = %e\n',j,d(j));
                end
            end
        end

        if (i_corr_fall ~= 2)
            alfa_d(j) = delta*alfa_d(j);
        end
        alfa = 0.0;
        return         
    end
end
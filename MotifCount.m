function [Label, Agr] = MotifCount(A,n,p,q,d)
    %See Section 9 of [GMPS19]
    B = A .* (A * A);
    phi_1 = asin(sqrt(betaincinv(p,d/2,1/2)));
    phi_2 = asin(sqrt(betaincinv(q,d/2,1/2)));
    ES = q*n + sqrt(6*q*n*log(n));
    ED = n*intersect(d+1,phi_1,phi_2)- sqrt(2*q*n*log(n));
    AA = (B >= ES | B <= ED);
    Label = assign(n,AA);
    Agr = abs([ones(1,n/2) -1*ones(1,n/2)] * Label / n);
    
    function label = assign(n,Adj)
        %Find the bipartition of a bipartite graph given adjacency matrix
        %If the graph is not bipartite, guarantee connected nodes are
        %assigned to the same cluster
        label = -1 * ones(n,1);
        visited = zeros(n,1);
        [~,ell] = max(sum(Adj));
        label(ell) = 1;
        visited(ell) = 1;
        neighbor = (Adj(:,ell) == 1);
        while min(visited - neighbor) < 0
            label(neighbor) = 1;
            visited(neighbor) = 1;
            neighbor = (max(Adj(:,neighbor),[],2) == 1);
        end
    end
        
    function area = intersect(Dim,t1,t2)
        %Find the surface 'area' of the intersection of two spherical caps
        if t1 >= 2*t2
            area = Cap(Dim,t2);
        elseif t1 + 2*t2 > 2*pi
            area = Cap(Dim,t1) - Cap(Dim,pi - t2);
        elseif (cos(t2))^2 >= cos(t1)
            tmin = atan(1/tan(t2) - cos(t1)/(cos(t2)*sin(t2)));
            area = Cap(Dim,t2) - J(Dim,tmin,t2) + J(Dim, t2+tmin, t1);
        else
            tmin = atan(-1/tan(t2) + cos(t1)/(cos(t2)*sin(t2)));
            area = J(Dim,t2-tmin,t1)+J(Dim,tmin,t2);
        end
        
        function areaa = Cap(D,t)
            %Find the surface 'area' of a spherical cap
            areaa = 0.5 * betainc((sin(t))^2,(D-1)/2,1/2);
        end
        
        function jay = J(D,tm,tM)
            %Numerical integration of a betainc function
            k = 10000;
            dt = (tM-tm)/k;
            t = tm+dt:dt:tM-dt;
            I = betainc(1-(tan(tm)./tan(t)).^2,(D-2)/2,1/2);
            J = sin(t).^(D-2).*I;
            J = sum(J*dt);
            jay = J * gamma(D/2) / (2 * sqrt(pi) * gamma((D-1)/2));
        end
    end
end
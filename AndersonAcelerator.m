function [ x, f_old, g_old, mAA, G, R, Q ] = AndersonAcelerator( mAA, g_curr, f_curr, iter, ...
                                             AAstart, G, f_old, g_old, beta, mMax, droptol, ...
                                             R, Q, minpvec, maxpvec )
%
    g_orig = g_curr;

    if AAstart == 1
        del_f = f_curr - f_old; del_g = g_curr - g_old;
        if mAA < mMax,
            G = [ G del_g ];
        else
            G = [ G(:,2:mAA) del_g ];
        end
        mAA = mAA + 1;
    end
    
    f_old = f_curr; g_old = g_curr;
    
    if mAA == 0
        % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
        x = g_curr;
    else
        % If mAA > 0, solve the least-squares problem and update the solution.
        if mAA == 1
            % If mAA == 1, form the initial QR decomposition.
            R(1,1) = norm(del_f);
            Q(:,1) = del_f/norm(del_f);
        else
            % If mAA > 1, update the QR decomposition.
            if mAA > mMax
                % If the column dimension of Q is mMax, delete the first column and
                % update the decomposition.
                [ Q, R ] = qrdelete(Q,R,1);
                mAA = mAA - 1;
                % The following treats the qrdelete quirk described below.
                % NÃO TEM ISSO NO ARTIGO!!!
                if size(R,1) ~= size(R,2),
                    Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                end
            end
            % Now update the QR decomposition to incorporate the new
            % column.
            for j = 1:mAA - 1
                R(j,mAA) = Q(:,j)'*del_f;
                del_f = del_f - R(j,mAA)*Q(:,j);
            end
            R(mAA,mAA) = norm(del_f);
            Q(:,mAA) = del_f/norm(del_f);
        end
        if droptol > 0
            % Drop residuals to improve conditioning if necessary.
            condDF = cond(R);
            while condDF > droptol && mAA > 1
                fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                [ Q, R ] = qrdelete(Q,R,1);
                G = G(:,2:mAA); % NÃO TEM ISSO NO ARTIGO!!!
                mAA = mAA - 1;
                % The following treats the qrdelete quirk described above.
                if size(R,1) ~= size(R,2),
                    Q = Q(:,1:mAA); R = R(1:mAA,:);
                end
                condDF = cond(R);
            end
        end
        % Solve the least-squares problem.
        gamma = R\(Q'*f_curr);
        % Update the approximate solution.
        x = g_curr - G*gamma;
        % Apply damping if beta is a function handle or if beta > 0
        % (and beta ~= 1).
        if isa(beta,'function_handle'),
            x = x - (1-beta(iter))*(f_curr - Q*R*gamma);
        else
            if beta > 0 && beta ~= 1 
                x = x - (1-beta)*(f_curr - Q*R*gamma);          
            end
        end
    end
    
    %---------------------------------------------------------------------%
    if min(g_old) < (min(minpvec) - 1e-9) || max(g_old) > (max(maxpvec) + 1e-9)
        posmin = find(g_old<(min(minpvec) - 1e-9));
        posmax = find(g_old>(max(maxpvec) + 1e-9));
        x(posmin) = g_orig(posmin);
        x(posmax) = g_orig(posmax);
    end
    %---------------------------------------------------------------------%
    
end


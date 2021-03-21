function modelupdate = find_stationary_adist(...
    p, model, income, grids, heterogeneity, varargin)
    % Finds the stationary distribution and transition matrix for a given
    % grids.a.vec.
    %
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    parser = inputParser;
    addParameter(parser, 'quiet', false);
    parse(parser, varargin{:});
    quiet = parser.Results.quiet;

    %% ----------------------------------------------------------------
    % FIND STATIONARY DISTRIBUTION
    % -----------------------------------------------------------------
    modelupdate = model;

    if ~quiet
        fprintf(' Computing state-to-state transition probabilities... \n');
    end

    nx = size(grids.a.vec, 1);
    R_bc = heterogeneity.R_broadcast;

    % Cash-on-hand as function of (a,yP,yF,yT)
    % start from generic 'a' distribution (even with returns het)
    x = grids.a.vec + income.netymat_broadcast;
    x = repmat(x, [1, 1, 1, p.nb]);
    
    % Saving interpolated onto this grid
    sav = zeros(nx,p.nyP,p.nyF,p.nb,p.nyT);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        x_iyP_iyF_ib = x(:,iyP,iyF,ib,:);
        sav_iyP_iyF_ib = model.savinterp{iyP,iyF,ib}(x_iyP_iyF_ib(:));
        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_ib, [nx 1 1 1 p.nyT]);
    end
    end
    end
    sav = max(sav, p.borrow_lim);

    % Transition matrix over (x,yP,yF,beta) full asset space
    modelupdate.statetrans = get_transition_matrix(p, income, grids,...
        nx, sav, R_bc);

    % Stationary distribution over states
    if ~quiet
        fprintf(' Finding ergodic distribution...\n');
    end
    q = get_distribution(p, grids, income, nx,...
        modelupdate.statetrans, heterogeneity, quiet);
%     [q,~] = eigs(modelupdate.statetrans',[],1,1);
%     q = q / sum(q(:));

    modelupdate.pmf = reshape(full(q'), [nx, p.nyP, p.nyF, p.nb]);

    tmp = reshape(full(q'), nx, []);
    modelupdate.pmf_a = sum(tmp, 2);
    
    % Get distribution over (x,yP,yF,beta)
    xdist = kron(income.yTdist, reshape(modelupdate.pmf, nx, []));
    modelupdate.xdist = reshape(xdist, [nx*p.nyT p.nyP p.nyF p.nb]);
    
    % Extend xvals to (nx*p.nyT,p.nyP,p.nyF,p.nyT)
    incvals = reshape(income.ymat, [p.nyP*p.nyF p.nyT]);
    incvals = permute(incvals, [2 1]);
    incvals = kron(incvals, ones(nx,1));
    incvals = reshape(incvals, [nx*p.nyT p.nyP p.nyF]);
    modelupdate.y_x = repmat(incvals, [1 1 1 p.nb]);
    modelupdate.nety_x = income.lumptransfer + (1-p.labtaxlow)*incvals...
        - p.labtaxhigh*max(incvals-income.labtaxthresh,0);
    modelupdate.nety_x = repmat(modelupdate.nety_x, [1 1 1 p.nb]);
    modelupdate.xvals = repmat(grids.a.vec, [p.nyT p.nyP p.nyF p.nb])...
        + modelupdate.nety_x;

    %% ----------------------------------------------------------------
    % POLICY FUNCTIONS ETC...
    % -----------------------------------------------------------------
    % Get saving policy function defined on xgrid
    modelupdate.sav_x = zeros(p.nx_DST*p.nyT,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP 
        modelupdate.sav_x(:,iyP,iyF,ib) = modelupdate.savinterp{iyP,iyF,ib}(...
            modelupdate.xvals(:,iyP,iyF,ib));
    end
    end
    end
    modelupdate.sav_x = max(modelupdate.sav_x, p.borrow_lim);

    % Policy functions associated with xdist
    savtax = p.compute_savtax(modelupdate.sav_x);
    modelupdate.con_x = modelupdate.xvals - modelupdate.sav_x ...
    	- savtax;
    
    % Mean assets
	modelupdate.mean_a = dot(modelupdate.pmf_a, grids.a.vec);

    if ~quiet
        fprintf(' A/Y = %2.5f\n', modelupdate.mean_a);
    end
end

%% ----------------------------------------------------------------
% TRANSITION MATRIX
% -----------------------------------------------------------------
function trans = get_transition_matrix(p, income, grids, nx, sav, R_bc)
	aprime_live = R_bc .* sav;

	% create interpolant object
    fspace = fundef({'spli', grids.a.vec, 0, 1});
    % get interpolated probabilities and take expectation over yT
    interp_live = 0;
    for k = 1:p.nyT
        ap_temp = aprime_live(:,:,:,:,k);
        interp_temp = reshape(funbas(fspace, ap_temp(:)), [], nx);
        interp_live = interp_live + income.yTdist(k) * interp_temp;
    end

    if p.Bequests
        interp_death = interp_live;
    else
        interp_death = sparse(nx*p.nyP*p.nyF*p.nb,nx);
        interp_death(:,grids.i0) = 1;
    end

    trans = sparse(nx*p.nyP*p.nyF*p.nb, nx*p.nyP*p.nyF*p.nb);
    col = 1;
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        transcol_live = kron(income.ytrans_live(:,col), ones(nx,1));
        transcol_death = kron(income.ytrans_death(:,col), ones(nx,1));
        
        transcol_live = transcol_live .* interp_live;
        transcol_death = transcol_death .* interp_death;

        % add new column to transition matrix
        trans(:,nx*(col-1)+1:nx*col) = ...
            (1-p.dieprob)*transcol_live + p.dieprob*transcol_death;
        col = col + 1;
    end
    end
    end
end

%% ----------------------------------------------------------------
% iTERATIVE METHOD TO FIND STATIONARY DISTRIBUTION
% -----------------------------------------------------------------
function q = get_distribution(p, grids, income, nx, statetrans,...
    heterogeneity, quiet)
	q = ones(nx,p.nyP,p.nyF,p.nb);

    % Create valid initial distribution
    yPdist = reshape(income.yPdist, [1,p.nyP]);
    yFdist = reshape(income.yFdist, [1,1,p.nyF,1]);
    zdist = reshape(heterogeneity.zdist, [1,1,1,heterogeneity.nz]);

    % Give households valid assets only (since with returns
    % heterogeneity and borrowing, lowest asset points cannot
    % be reached or left by some households)
    if numel(p.r) > 1
        for ib = 1:p.nb
            tmp1 = sum(reshape(q(:,:,:,ib), [], 1));
            for ia = 1:p.nx
                if grids.a.vec(ia) < grids.s.vec(1) * p.R(ib)
                    q(ia,:,:,ib) = 0;
                else
                    break
                end
            end
            tmp2 = sum(reshape(q(:,:,:,ib), [], 1));
            q(:,:,:,ib) = q(:,:,:,ib) * tmp1 / tmp2;
        end
    end
    q = q .* yPdist .* yFdist .* zdist;
    q = q(:)' / sum(q(:));

    diff = 1; 
    iter = 1;
    while diff>1e-9 && iter < 5e5
        z = q * statetrans;
        diff = norm(z-q);
        q = z;
        
        if ~quiet && (mod(iter,500) == 0)
            fprintf('  Diff = %5.3E, Iteration = %u \n',diff,iter);
        end
        iter = iter + 1;
    end

    if iter >= 5e5
        error('No conv to statdist, diff = %5.3e',diff)
    end
end
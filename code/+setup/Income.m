classdef Income < handle
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    properties (SetAccess = private)
        % Read the persistent income process from mat file, boolean
        Load_yP;

        % Read the transitory income process from mat file, boolean
        Load_yT;

        % Structure containing data read from mat file
        ImportedVariables;
        
        % Parameters object
        p;

        % Transition matrix of z-dimension (e.g. discount factor het)
        ztrans;
        
        % Persistent income component
        logyPgrid;
        yPdist;
        yPtrans;
        yPgrid
        yPcumdist;
        yPcumtrans;
        
        % Transitory income component
        logyTgrid;
        yTdist;
        yTgrid;
        yTcumdist;
        
        % Fixed income component
        logyFgrid;
        yFgrid;
        yFdist;
        yFcumdist;

        % Transition matrix over (yF, yP) states
        ytrans;

        % One-period mean of gross labor income
        meany1;

        % Matrix of all possible gross income values
        ymat;

        % Distribution over all possible income values
        ymatdist;

        % Threshold above which labor is subject to tax
        labtaxthresh;

        % One-period value of lump sum transfer
        lumptransfer;
        
        % Net income values
        netymat;
        netymat_broadcast;
        netymatEGP;
        netymatDST;
        meannety1;
        minnety;
        
        % Share of households for which net labor income exceeds gross labor income
        fraction_net_transfer;
        
        % Transition matrix across income and z-dimension states, conditional on life/death
        ytrans_live;
        ytrans_death;
    end
    
    methods
        function obj = Income(p, heterogeneity)
            obj.p = p;
            obj.ztrans = heterogeneity.ztrans;

            % Read persistent process if an income path is specified
            obj.Load_yP = ~isempty(p.IncomeProcess);
            obj.Load_yT = false;
            if obj.Load_yP
                obj.ImportedVariables = load(p.IncomeProcess);
                if ~isempty(obj.ImportedVariables.logyTgrid)
                    obj.Load_yT = true;
                end
            end
            
            % Computations
            obj.get_persistent_income();
            obj.get_transitory_income();
            obj.get_fixed_effect();
            obj.get_other_income_variables();
            
            % Error checks
            if size(obj.yTgrid,2)>1 || size(obj.yFgrid,2)>1 || size(obj.yPgrid,2)>1
                error('All income grids must be column vectors')
            end
            if size(obj.yTdist,2)>1 || size(obj.yFdist,2)>1 || size(obj.yPdist,2)>1
                error('All income distributions must be column vectors')
            end
        end
        
        function get_persistent_income(obj)
            if obj.Load_yP
                obj.logyPgrid = obj.ImportedVariables.logyPgrid;
                obj.yPdist = obj.ImportedVariables.yPdist;
                obj.yPtrans = obj.ImportedVariables.yPtrans;
                obj.p.nyP = length(obj.logyPgrid);
                obj.logyPgrid = reshape(obj.logyPgrid,[],1);
                obj.yPdist = reshape(obj.yPdist,[],1);
            elseif obj.p.nyP > 1
                [obj.logyPgrid, obj.yPtrans, obj.yPdist] ...
                    = aux.rouwenhorst(obj.p.nyP, -0.5*obj.p.sd_logyP^2, obj.p.sd_logyP, obj.p.rho_logyP);
            else
                obj.logyPgrid = 0;
                obj.yPdist = 1;
                obj.yPtrans = 1;
            end
            
            obj.yPgrid = exp(obj.logyPgrid);
            obj.yPgrid = obj.yPgrid/(obj.yPdist'*obj.yPgrid);
            obj.logyPgrid = log(obj.yPgrid);
            obj.yPcumdist = cumsum(obj.yPdist,1);
            obj.yPcumtrans = cumsum(obj.yPtrans,2);
        end
        
        function get_transitory_income(obj)
            if obj.Load_yT
                obj.logyTgrid = obj.ImportedVariables.logyTgrid;
                obj.yTdist = obj.ImportedVariables.yTdist;
                obj.p.nyT = length(obj.logyTgrid);
                obj.logyTgrid = reshape(obj.logyTgrid,[],1);
                obj.yTdist = reshape(obj.yTdist,[],1);
                obj.yTcumdist = cumsum(obj.yTdist,1);
            elseif (obj.p.nyT>1) && (obj.p.sd_logyT>0)
                %moments of mixture distribution
                lmu2 = obj.p.lambdaT.*obj.p.sd_logyT^2;
                lmu4 = 3.*obj.p.lambdaT.*(obj.p.sd_logyT^4);

                %fit those moments
                mu1 = 0;
                optionsNLLS = optimoptions('lsqnonlin','Display','Off');
                lpar = lsqnonlin(@(lp) aux.discretize_normal_var_kurt(...
                    lp,obj.p.nyT,mu1,lmu2,lmu4,obj.p.sd_logyT),[0.1 1],[],[],optionsNLLS);
                [lf,lx,lp] = aux.discretize_normal_var_kurt(lpar,obj.p.nyT,mu1,lmu2,lmu4,obj.p.sd_logyT);
                obj.logyTgrid = lx;
                obj.yTdist = lp;
                obj.yTcumdist = cumsum(obj.yTdist,1);

            elseif obj.p.nyT==1
                obj.logyTgrid = 0;
                obj.yTdist = 1;
                obj.yTcumdist = 1;
            elseif obj.p.sd_logyT==0
                obj.logyTgrid = linspace(-2, 2, obj.p.nyT)';
                obj.yTdist = zeros(obj.p.nyT, 1);
                obj.yTdist(median(1:obj.p.nyT)) = 1;
                obj.yTcumdist = cumsum(obj.yTdist);
            end
            
            obj.yTgrid = exp(obj.logyTgrid);
            obj.yTgrid = obj.yTgrid / (obj.yTdist'*obj.yTgrid);
            obj.logyTgrid = log(obj.yTgrid);
        end
        
        function get_fixed_effect(obj)
            if obj.p.nyF>1
                width = fzero(@(x) aux.discrete_normal(obj.p.nyF,-0.5*obj.p.sd_logyF^2 ,obj.p.sd_logyF ,x),2);
                [~,obj.logyFgrid,obj.yFdist] = aux.discrete_normal(obj.p.nyF,-0.5*obj.p.sd_logyF^2 ,obj.p.sd_logyF ,width);
                obj.logyFgrid = reshape(obj.logyFgrid,[],1);
                obj.yFdist = reshape(obj.yFdist,[],1);
            elseif obj.p.nyF==1
                obj.logyFgrid = 0;
                obj.yFdist = 1;
            end
            obj.yFgrid = exp(obj.logyFgrid);

            % normalize fixed effect such that mean gross y = 1 if annual, 1/4 if quarterly
            obj.yFgrid = obj.yFgrid/(obj.yFdist'*obj.yFgrid*obj.p.freq);
            obj.logyFgrid = log(obj.yFgrid);
            obj.yFcumdist = cumsum(obj.yFdist,1);
        end
        
        function get_other_income_variables(obj)
            % transition probabilities for yP-yF combined grid
            obj.ytrans = kron(eye(obj.p.nyF), obj.yPtrans);

            % construct matrix of y combinations
            obj.ymat = repmat(obj.yPgrid,obj.p.nyF,1) ...
            	.* kron(obj.yFgrid,ones(obj.p.nyP,1)) * obj.yTgrid';

            % distribution of ymat
            obj.ymatdist = repmat(obj.yPdist,obj.p.nyF,1) ...
            	.* kron(obj.yFdist,ones(obj.p.nyP,1)) * obj.yTdist';

            % find mean y
            % isolate unique (yT,yF,yP) combinations
            temp = sortrows([obj.ymat(:) obj.ymatdist(:)],1);
            ysort = temp(:,1);
            ysortdist = temp(:,2);
            ycumdist_sort = cumsum(ysortdist);

            % 1-period statistics
            obj.meany1 = obj.ymat(:)' * obj.ymatdist(:);
            totgrossy1 = obj.meany1;

            % find tax threshold on labor income
            if numel(ysort)>1
                obj.labtaxthresh = lininterp1(ycumdist_sort,ysort,obj.p.labtaxthreshpc);
            else
                obj.labtaxthresh = 0;
            end    

            % find net income
            totgrossyhigh = max(obj.ymat(:)-obj.labtaxthresh,0)' * obj.ymatdist(:);
            % obj.lumptransfer = obj.p.labtaxlow * totgrossy1 ...
            % 	+ obj.p.labtaxhigh * totgrossyhigh;
            obj.lumptransfer = obj.p.lumptransfer;

            % netymat is N by nyT matrix
            obj.netymat = obj.lumptransfer + (1-obj.p.labtaxlow) * obj.ymat ...
            	- obj.p.labtaxhigh * max(obj.ymat-obj.labtaxthresh,0);
            obj.meannety1 = obj.netymat(:)' * obj.ymatdist(:);

            obj.minnety = min(obj.netymat(:));
            
            % fraction of households that recieve net transfer from gov
            obj.fraction_net_transfer = (obj.netymat(:)>obj.ymat(:))' * obj.ymatdist(:);

            % net y values on HJB grid
            obj.netymat_broadcast = reshape(obj.netymat,[1 obj.p.nyP obj.p.nyF 1 obj.p.nyT]);
            obj.netymatEGP = repmat(obj.netymat_broadcast,[obj.p.nx 1 1 obj.p.nb 1]);
            obj.netymatDST = repmat(obj.netymat_broadcast,[obj.p.nx_DST 1 1 obj.p.nb 1]);

            % full transition matrix with beta and IES transitions, excluding and including death
            if obj.p.ResetIncomeUponDeath
                yPtrans_death = repmat(obj.yPdist',obj.p.nyP,1);
            else
                yPtrans_death = obj.yPtrans;
            end

            obj.ytrans_live = kron(obj.ztrans,kron(eye(obj.p.nyF),obj.yPtrans));
            obj.ytrans_death = kron(obj.ztrans,kron(eye(obj.p.nyF),yPtrans_death));
        end
    end
end
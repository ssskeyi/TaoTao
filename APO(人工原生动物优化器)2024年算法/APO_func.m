%%
% Artificial Protozoa Optimizer (APO): A novel bio-inspired metaheuristic algorithm for engineering optimization
function [bestProtozoa,bestFit,cuver] = APO_func(pop_size,iter_max,Xmin,Xmax,dim,fobj)
% random seeds
stm = RandStream('swb2712','Seed',sum(100*clock));
RandStream.setGlobalStream(stm);
% global best
Xmin = Xmin.*ones(1,dim);
Xmax = Xmax.*ones(1,dim);
ps = pop_size;  % ps denotes protozoa size
np = 1;         % np denotes neighbor pairs     np_max can be set in floor((ps-1)/2)
pf_max = 0.1;   % pf_max denotes proportion fraction maximum
%
tic;
protozoa=zeros(ps,dim);    % protozoa
newprotozoa=zeros(ps,dim); % new protozoa
epn=zeros(np,dim); % epn denotes effect of paired neighbors
% initilization
for i = 1:ps
    protozoa(i,:) = Xmin + rand(1,dim).*(Xmax-Xmin);
    protozoa_Fit(i) = fobj(protozoa(i,:));
end

% find the bestProtozoa and bestFit
[bestval,bestid] = min(protozoa_Fit);
bestProtozoa = protozoa(bestid,:);  % bestProtozoa
bestFit = bestval; % bestFit
cuver(1) = bestFit;
%%  Main loop
for iter=2:iter_max
    [protozoa_Fit,index] = sort(protozoa_Fit);
    protozoa= protozoa(index,:);
    pf = pf_max*rand; % proportion fraction
    ri=randperm(ps,ceil(ps*pf)); % rank index of protozoa in dormancy or reproduction forms
    for i=1:ps
        if ismember(i,ri) %  protozoa is in dormancy or reproduction form
            pdr=1/2*(1+cos((1-i/ps)*pi)); % probability of dormancy and reproduction
            if rand<pdr  % dormancy form
                newprotozoa(i,:)=  Xmin + rand(1,dim).*(Xmax-Xmin);
            else  % reproduction form
                flag=[1,-1];  % +- (plus minus)
                Flag=flag(ceil(2*rand));
                Mr=zeros(1,dim); % Mr is a mapping vector in reproduction
                Mr(1,randperm(dim,ceil(rand*dim)))=1;
                newprotozoa(i,:)= protozoa(i,:) + Flag*rand*(Xmin+rand(1,dim).*(Xmax-Xmin)).*Mr;
            end
        else  % protozoa is foraging form
            f= rand*(1+cos(iter/iter_max*pi)); % foraging factor
            Mf=zeros(1,dim);  % Mf is a mapping vector in foraging
            Mf(1,randperm(dim,ceil(dim*i/ps)))=1;
            pah= 1/2*(1+cos(iter/iter_max*pi)); % probability of autotroph and heterotroph
            if rand<pah  % protozoa is in autotroph form
                j= randperm(ps,1); % j denotes the jth randomly selected protozoa
                for k=1:np % np denotes neighbor pairs
                    if i==1
                        km=i; % km denotes the k- (k minus)
                        kp=i+randperm(ps-i,1); % kp denotes the k+ (k plus)
                    elseif i==ps
                        km=randperm(ps-1,1);
                        kp=i;
                    else
                        km=randperm(i-1,1);
                        kp=i+randperm(ps-i,1);
                    end
                    % wa denotes weight factor in the autotroph forms
                    wa=exp(-abs(protozoa_Fit(1,km)/(protozoa_Fit(1,kp)+eps)));
                    % epn denotes effect of paired neighbors
                    epn(k,:)=wa*(protozoa(km,:)-protozoa(kp,:));
                end
                newprotozoa(i,:)= protozoa(i,:)+ f*(protozoa(j,:)-protozoa(i,:)+1/np*sum(epn,1)).*Mf;

            else   % protozoa is in heterotroph form
                for k=1:np % np denotes neighbor pairs
                    if i==1
                        imk=i;   % imk denotes i-k (i minus k)
                        ipk=i+k; % ipk denotes i+k (i plus k)
                    elseif i==ps
                        imk=ps-k;
                        ipk =i;
                    else
                        imk=i-k;
                        ipk=i+k;
                    end
                    % neighbor limit range in [1,ps]
                    if  imk<1
                        imk=1;
                    elseif ipk>ps
                        ipk=ps;
                    end
                    % denotes weight factor in the heterotroph form
                    wh=exp(-abs(protozoa_Fit(1,imk)/(protozoa_Fit(1,ipk)+eps)));
                    epn(k,:)=wh*(protozoa(imk,:)-protozoa(ipk,:));
                end
                flag=[1,-1];  % +- (plus minus)
                Flag=flag(ceil(2*rand));
                Xnear=(1+Flag*rand(1,dim)*(1-iter/iter_max)).* protozoa(i,:);
                newprotozoa(i,:)=protozoa(i,:)+f*(Xnear-protozoa(i,:)+1/np*sum(epn,1)).*Mf;
            end
        end
        newprotozoa(i, : ) = Bounds( newprotozoa(i, : ), Xmin, Xmax );
        newprotozoa_Fit(i)= fobj(newprotozoa(i, : ));

    end

    bin = (protozoa_Fit > newprotozoa_Fit)';
    protozoa(bin==1,:) = newprotozoa(bin==1,:);
    protozoa_Fit(bin==1) = newprotozoa_Fit(bin==1);
    [bestFit,bestid] = min(protozoa_Fit);
    bestProtozoa = protozoa(bestid,:);
    cuver(iter)=bestFit;
end
recordtime = toc;
end

% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
% Apply the lower bound vector
temp = s;
I = temp < Lb;
temp(I) = Lb(I);

% Apply the upper bound vector
J = temp > Ub;
temp(J) = Ub(J);
% Update this new move
s = temp;
end


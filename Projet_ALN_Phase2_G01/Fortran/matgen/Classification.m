
%Number of simulations
Nens=20;
%Quality of the subspace approximation
percentInfo=0.95;

%Reading the observations
%%%
% TODO: select the file associated to your group number  
load('observation.mat');
	  
% Initiatlization	  
Pi=zeros(3,1);

for GWi = 1:3
    %Generation of the data set
    Fi = Model(GWi,Nens);
    
    %Computation of the mean and anomalies 
    mFi= mean(Fi,2);
    Zi=Fi-repmat(mFi,1,Nens);
    
    %%%
    [Ui,Si,Vi] = svd(Zi,0);
    DS = diag(Si);
    if (DS(1)==0)
      disp('Matrix null')
      return
    end

    %%%%
    % Select the vectors associated with the most dominant singular values.
    % This is done accordingly to Equation (1).
    normZ2=norm(Zi,'fro')^2;
    converged=1;
    n=length(Zi(:,1));
    while ((DS(converged)/DS(1)>1-percentInfo)&&(converged<n+1)) 
      converged=converged+1;
    end
    converged=converged-1; 
       
    Ui = Ui(:,1:converged);
    Vi = Vi(:,1:converged);
    %%%

    Zobs=Fobs-mFi;
    tmp=(Ui')*Zobs; % to prevent "out of memory" issues 
    Pi(GWi)=norm(Zobs-Ui*tmp);
    
end

figure(1)
bar(Pi)

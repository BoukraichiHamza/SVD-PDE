%%% Set check to true for validation %%%
check=false;

% Number of simulations
%%%
% TODO : assess the impact fot this parameter on the results
% of the reconstruction
Nens = 50;

% Wind parameter
GW = 1;
% Stopping criterions 
maxIter=100;
epsilon=1.e-8;
% TODO : assess the impact of this parameter on the results
% of the reconstruction
percentInfo = 0.95;

% Generate the simulations
F = Model(GW,Nens);

% Ensemble mean
muF = mean(F,2);
% Compute the anomaly matrix
Z   = F - repmat(muF,1,Nens);

%%%%%%%  Compute the SVD of A    %%%%%%%
if (check)   
  tic;
  [U,S,V] = svd(Z,0);
  d = diag(S);
  if (d(1)==0)
    disp('Alert: the matrix is null')
    return
  end
 
  % Trace of Z*Z'
  n=length(Z(:,1));  
  converged=1;
  while ((d(converged)/d(1)>1-percentInfo)&&(converged<n+1)) 
    converged=converged+1;
  end
  converged=converged-1; 
      
  U = U(:,1:converged);
  V = V(:,1:converged);
  time=toc;
  fprintf(['%d singular values were found in %7.3f seconds:\n'],converged,time);

else

  % Initialization: some subspace iteration method parameters.
  % Rem: these initial values are not.meant to be relevant. 
  p=2;
  m=45;
  
  % Timer
  tic;
  
  % Trace of Z*Z'
  normZ2=norm(Z,'fro')^2;
  
  n=length(Z(:,1));  
  Vs=zeros(n,m);
  for k=1:m
    Vs(k,k)=1;
  end
  s=zeros(m,1);
  
  converged=0;
  iter=0;  
  condition=0;
  
  while((converged<m)&&(iter<maxIter)) 
    iter=iter+1;
    
    %Compute Y=(Z*Z')^p*V
    Y=Vs;
    for k=1:p
      Y=Z*((Z')*Y);
    end
   
    %Gram-Schmidt orthogonalization
    Vs=Y;
    for i=1:m
       for j=1:i-1
    	  tmp=Vs(:,j)'*Vs(:,i);
          Vs(:,i) = Vs(:,i) - tmp*Vs(:,j);
       end
       Vs(:,i) = Vs(:,i)/norm(Vs(:,i));
    end
  
    %Rayleigh quotient
    H=(Vs')*(Z*((Z')*Vs));
  
    %eig-decomposition of H
   [Ens,dns]=eig(H);
   [d,index]=sort(diag(dns),'descend');
    E=Ens(:,index);
   
    %V=V*X
    Vs=Vs*E; 
  
    % Check wich are the eigenvalues that have converged 
    conv=0;
    for i=converged+1:m
      res=norm(Z*((Z')*Vs(:,i))-d(i)*Vs(:,i),'fro')/normZ2; 
      if (res>epsilon)
        break
      end
      conv=conv+1;
      s(i)=sqrt(d(i));  
      condition=1-s(i)/s(1);  
    end
    converged=converged+conv;
  
    %Stopping criterion
    if(condition>=percentInfo)
      disp('Convergence')
      converged=max(converged-1,1);
      break
    end  
    
  end 
 
  %Singular vectors
  if(converged==0)
     disp('maxIter reached without convergence')
     return
  else   
     U=Vs(:,1:converged);
     V=Z'*U;
     for k=1:converged
        V(:,k)=V(:,k)/s(k);
     end
  end
  
  %End timer 
  time=toc;
  
  fprintf(['%d singular values were found in %7.3f seconds.\n'],converged,time);
end

%%%%%%%       Reconstruction        %%%%%%%
[X, ns, nt] = Model(GW,1);
X0 = X(1:ns,:);
%%%%
% Reconstruct X with X0
Zp = zeros(size(X)); 
%%%%
Z0=X0-muF(1:ns);
alpha=(U(1:ns,:)'*U(1:ns,:))\(U(1:ns,:)'*Z0);
Zp=U*alpha;

%%%% Compute the error %%%%
Xp = Zp + muF;
error=norm(Xp-X)/norm(X);
fprintf('error = %f\n',error);

%%%% Display %%%%
global Lx Ly Nx Ny;

% draw result
    x = linspace(0,Lx,Nx);     %  Independent variable x
    y = linspace(0,Ly,Ny);     %  Independent variable y
    [Mx, My] = meshgrid(x,y);  %  2D arrays, mainly for plotting
    Mx = Mx'; My = My';        %  MatLab is strange!

      figure(2)
     
 for tt=1:nt
          set(gcf,'Renderer','Painters')
          subplot(1,2,1);
          z = X((tt-1)*ns+1:tt*ns,1);
          z = reshape(z,Nx,Ny);
          surf(Mx,My,z); shading('interp');
          axis([0,Lx,0,Ly ,5000,6000]);
	  pbaspect([3 1 3])
	  title('Solution')
	  
      
          subplot(1,2,2);
          zappr = Xp((tt-1)*ns+1:tt*ns,1);
          zappr = reshape(zappr,Nx,Ny);
          h2=surf(Mx,My,zappr); shading('interp');
	  axis([0,Lx,0,Ly ,5000,6000]);
	  pbaspect([3 1 3])
	  title('Reconstruction') 
         drawnow	 
  end

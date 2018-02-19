%%% Set check to true for validation %%%
check=false;

% Number of simulations
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
  p=2; % Block iterations
  m=45;  % Maximum dimension of the subspace 
  
  % Timer
  tic;
  
  converged=0; % Number of singular values that have been found
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% TO DO : computation of the subspaces based on the power iteration method
  
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
    
    if mod(iter,2) == 0
        %Compute Y=(Z*Z')^p*V
        Y1=Vs;
        for k=1:p
            if mod(k,2)
                Y1=(Z')*Y1;
            else
                Y1=Z*Y1;
            end
        end
        
    else 
        %Compute Y=(Z'*Z)^p*U
        Y2=Us;
        for k=1:p
        if mod(k,2)
                Y2=Z*Y2;
            else
                Y2=(Z')*Y2;
            end
        end
    end
    
    if mod(p,2) == 1
       if mod(iter,2) == 0
           Y1 = Z*Y1;
       else 
           Y2 = Z'*Y2;
       end
    end
    
    %Gram-Schmidt orthogonalization
     if mod(iter,2) == 0
           M = Y1;
       else 
           M = Y2;
     end
     
    for i=1:m
       for j=1:i-1
    	  tmp=M(:,j)'*M(:,i);
          M(:,i) = M(:,i) - tmp*M(:,j);
       end
       M(:,i) = M(:,i)/norm(M(:,i));
    end
  
    %Rayleigh quotient
    if mod(iter,2) == 0
           H=(M')*(Z*((Z')*M));
       else 
           H=(M')*((Z')*(Z*M));
     end
    
  
    %eig-decomposition of H
   [Ens,dns]=eig(H);
   [d,index]=sort(diag(dns),'descend');s=zeros(m,1);
    E=Ens(:,index);
   
    %M=M*X
    M=M*E;
    
     if mod(iter,2) == 0
         %Save V and compute U = ZV
         Vs = M;
         Us = Z'*M;
        
    else 
        %save U and compute V = Z'U
        Us = M;
        Vs = Z*M;
    end
    
  
    % Check wich are the eigenvalues that have converged 
    conv=0;
    for i=converged+1:m
        if mod(iter,2) == 0
           res=norm(Z*((Z')*M(:,i))-d(i)*M(:,i),'fro')/normZ2; 
       else 
          res=norm((Z')*(Z*M(:,i))-d(i)*M(:,i),'fro')/normZ2; 
     end
      
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
    iter=iter+1;
  end 
  
   %Singular vectors
  if(converged==0)
     disp('maxIter reached without convergence')
     return
  else   
     U=Us(:,1:converged);
     V=Vs(:,1:converged);
  end
  
  % END TO DO
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
alpha=(V(1:ns,:)'*V(1:ns,:))\(V(1:ns,:)'*Z0);
Zp=V*alpha;

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

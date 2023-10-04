function [patameters_PLS]=pls_svd(M,Y,Fac)
%%%  
%%%
%%%
%%% Partial Least Squares: SVD
%%%  
%%%
%%%
%     Inputs:  
%     M    : block A matrix (number of samples  x dim1) - zero mean !
%     Y    : block B matrix (number of samples  x dim2) - zero mean !
%     Fac  : number of latent vectors (components)  to extract 

%     Outputs: 
%     T    : matrix of latent vectors (number of samples x Fac)  
%     V    : weight maxtrix           
        
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
  
[~,cX]=size(M);
[~,cY]=size(Y);
T=zeros(cX,Fac);
V=zeros(cX,Fac);
C=zeros(cY,Fac);
Q=zeros(cY,Fac);

fprintf('Train the PLS-based TVS with dimension %d:\n',Fac);
for i = 1 :Fac
%     fprintf('%i ',i);
    MY=M'*Y;
    [v,~,c]=svds(MY,1);
    w=M*v;
    u=Y*c;
    t=M'*w/(w'*w);
    q=Y'*u/(u'*u);
    b=u'*w/(w'*w);
    V(:,i)=v;
    C(:,i)=c;
    T(:,i)=t;
    Q(:,i)=q;
    B(i)=b;
    %fprintf('Loading matrix get\n');
    M = M - w*t';
    Y = Y - b*w*q';
    %fprintf('Residual matrix get\n');
end

%     XY=X'*Y;
%     [W,S,C]=svds(XY,Fac);
%     for i = 1 :Fac
%         t=X*W(:,i);
%         u=Y*C(:,i);
%         p=X'*t/(t'*t);
%         q=Y'*u/(u'*u);
%         b=u'*t/(t'*t);
% %         T(:,i)=t;
%         P(:,i)=p;
%     end

patameters_PLS.T = T;
patameters_PLS.V = V;
patameters_PLS.Q = Q;
patameters_PLS.C = C;
patameters_PLS.B = B;




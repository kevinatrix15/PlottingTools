function[F, Nx, Ny, Nxy]=LaplacianIntegrator(bedwidth_x,bedwidth_y,x_0,y_0,Nodes)
% clear all;
% clc;
format long
syms x y;
% Nodes is the node number we are going to choose to make a standard mesh
% for this calculation
%Nodes=8;

%bedwidth_x=1e-3;
%bedwidth_y=1e-3;

state_res=1;

h(1)=(1/state_res)*bedwidth_x/(Nodes-1);
h(2)=(1/state_res)*bedwidth_y/(Nodes-1);
np(1)=state_res*(Nodes-1)+1;
np(2)=state_res*(Nodes-1)+1;


[cord_state,elem_state]=box_creater(np(1),np(2),h(1),h(2));

Perm_mat=zeros(Nodes+2,2*Nodes+2);
Perm_mat(1,1)=1;
Perm_mat(Nodes+2,2*Nodes+2)=1;
k=2;
for i=2:Nodes+1
    Perm_mat(i,k)=1;
    Perm_mat(i,k+1)=1;
    k=k+2;
end

Lxfull=(1/h(1))*Perm_mat*...
  kron(...
    eye(Nodes+1,Nodes+1),[1 -1;-1 1]...
  )*Perm_mat';
Lx=Lxfull(2:Nodes+1,2:Nodes+1);

Myfull=(1/6)*h(1)* Perm_mat*...
  kron(...
    eye(Nodes+1,Nodes+1),[2 1;1 2]...
  )*Perm_mat';
My=Myfull(2:Nodes+1,2:Nodes+1);

[V,D]=eig(Lx,My);

V_kron = kron(V,V);
Vt_kron = kron(V', V');


Diag=kron(eye(Nodes,Nodes),D)+kron(D,eye(Nodes,Nodes));
Diag_sq=Diag^2;

M_isp=(My^-1)*(My^-1);
M_isp_sq= kron(M_isp,M_isp);

%x_0,y_0 read from main code tells the laser location 
g= @(x,y) 4*exp(-(x-x_0)^2-(y-y_0)^2)*(x^2-2*x*x_0 + x_0^2 + y^2 -2*y*y_0 + y_0^2 - 1);
% gg=zeros(Nodes^2,1);
% gg(Nodes^2,1)=1;
gg=eval(subs(g,{x,y},{cord_state(:,1),cord_state(:,2)}));
F=-M_isp_sq*(Diag_sq^(-1)*gg);

Nx = -kron(Lx, eye(Nodes,Nodes))*F/(h(2)^2);
Ny = -kron(eye(Nodes,Nodes), Lx)*F/(h(1)^2);

A = zeros(Nodes,Nodes);
A(1,1) = 1;
k = 1;
for i = 2:Nodes
    A(i,k) = -1;
    A(i,k+1) = 1;
    k = k+1;
end

    
Nxy = kron(A, eye(Nodes,Nodes))*kron(eye(Nodes,Nodes),A')*F/(h(1)*h(2))



clc
clear all
syms GJ_root GJ_tip GJ b y w(y)  t(Y) alphae(y) M(y) Pdyn l(y) P(y) c(y) w(y) theta(y) b1 b2 b3 b4 b5 a1 a2 a3 a4 a5 y z

[Mc,Cla,Clac,Cmac,sweep,s,N,g,W,elevation,rho,T]=deal(0.8,6.28,5.05,-0.015,34*pi/180,21.33,3,9.81,66000,6000,0.54,251.7);
%Pdyn=(sqrt(1.4*287*T)*Mc)^2/2*.66; %N/m^2
%LINEARLY VARYING TERMS
[EI(y),GJ(y),c(y),M(y),ec(y),x_cg(y)]=deal(4e8-.1782e8*y,2.5e8-.109e8*y,5.285-.144253*y,800-28.13*y,.68705-0.018753*y,-0.106+.00288*y);

% no of terms
[n, m, w(y), theta(y)] = deal(5, 5, 0,0);
[a,b]=deal(sym('a',[1 n]),sym('b',[1 m]));
for i=1:n
w(y)=w(y)+a(i)*(y^4+6*s^2*y^2-4*y^3*s)^i;
end
for i=1:m
theta(y)=theta(y)+ b(i)*(-2*s*y+y^2)^i;
end
syms alpha
alphae(y)= alpha + theta(y)*cos(sweep)-diff(w(y),y)*sin(sweep);
syms l(y) m(y)
l(y)= Pdyn*Clac*c(y)*alphae(y);
P(y)= (l-N*M(y)*g)*cos(sweep);
m(y)=(l(y)*ec(y)+Pdyn*((c(y))^2)*Cmac-N*M(y)*g*x_cg(y))*(cos(sweep))^2;
D=EI(y)*diff(diff(w(y),y),y);
%residuals
R1(y)=diff(diff(D,y),y)-P(y)-diff((l(y)*ec(y)+Pdyn*((c(y))^2)*Cmac-N*M(y)*g*x_cg(y)),y)*sin(sweep)*cos(sweep);
R2(y)=diff((GJ(y)*diff(theta(y),y)),y)+m(y);
[Q1,Q2] = deal(R1(y)*R1(y),R2(y)*R2(y));
% using weighted integral residual -least square error
[I1,I2] = deal(int(Q1, 0, s),int(Q2, 0, s));
eq = cell(10,1);
eq{1}=diff(I1 , a1);
eq{2}=diff(I1 , a2);
eq{3}=diff(I1 , a3);
eq{4}=diff(I1 , a4);
eq{5}=diff(I2 , b1);
eq{6}=diff(I2 , b2);
eq{7}=diff(I2 , b3);
eq{8}=diff(I2 , b4);
eq{9}=diff(I1 , a5);
eq{10}=diff(I2 , b5);
% divergence matrix formation 
a= [a1;a2;a3;a4;b1;b2;b3;b4;a5;b5];
matrix = cell(10,1);
for i=1:10
    matrix{i} = fliplr(coeffs(eq{i},a));
end
%calculation of rigid angle of attack
L = vpa(2*int(l(y),[0,s]),7);
alpha_val = subs(L-N*W*g);
alpha = vpa(solve(alpha_val),7);
D= subs(matrix);
coeff = D(:,1:10); 
% Solving for divergence dynamic pressure
syms Pdyn
eqn = det(coeff) == 0;
PdynD = solve(eqn, Pdyn);
% Converting symbolic expression to decimal 
PdynD = vpa(subs(PdynD), 5);
% Removing imaginary and negative values
PdynD = PdynD((imag(PdynD)==0) & (PdynD>=0));
% minimum positive value for divergence pressure
PdynD = min(PdynD);
% Calculating divergence velocity and Mach number
Divergence_speed = sqrt(2*PdynD/rho^2);
MachNo = Divergence_speed/318; % Divergence Mach no at 6000m elevation
% Print results
fprintf('\n Divergence dynamic pressure = %0.4f Pa', PdynD);
fprintf('\n Divergence Velocity = %0.4f m/s', Divergence_speed);
fprintf('\n Divergence Mach No. = %0.4f', MachNo);

%-----------------------------------------------------------------------------------%


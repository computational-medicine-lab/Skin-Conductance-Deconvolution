function [A, B] = create_A_B_matrix_ss_multires(tau, Nu, Fsu, Fsy)

% theta, skin
% Nu
% Fsu
% Fsy

Fs = Fsu;
Ts = 1/Fs;
A = [-1/tau(1) 0; +1/tau(2) -1/tau(2)];      B = [1/tau(1) ;0];  C = [0 1];  D = 0;
 sys = ss(A,B,C,D); sysd = c2d(sys, Ts);
 Ad = sysd.A; Bd = sysd.B; Cd = sysd.C; Dd = sysd.D; 
 
 y0 = 0.5; r0 = 0; x = [];
 x0 = [r0; y0];
 x(:,1) = x0;
 y(1) = Cd * x(:,1);
 
 D1 = zeros(Nu,Nu);

 temp = 1;
 for i=1:Nu
     F1(i,:) = Cd * Ad^(i-1);
 end
 temp = 1;
 for i=1:Nu
     D_(Nu-i+1) = Cd * temp * Bd;
     temp = temp * Ad;
 end
 for i=1:Nu
     D1(Nu-i+1,:) =  [D_(i+1:end) zeros(1,i)];
 end
 
 downsampling_factor = Fsu/Fsy;
 
 A = downsample(F1,downsampling_factor);
 B = downsample(D1,downsampling_factor);
 
end
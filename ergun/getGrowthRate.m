function r = phi(n)
   porosity = 0.8;
   k = sqrt(2*4*pi^2);
   rho1 = 1000;
   rho2 = 1000*porosity+8000*(1-porosity);
   mu1 = 0.218375*rho1;
   mu2 = mu1*(1-(1-porosity)/0.63)^(-2.5*0.63);
   alpha1 = rho1/(rho1+rho2);
   alpha2 = rho2/(rho1+rho2);
   v1 = mu1/rho1;
   v2 = mu2/rho2;
   q1 = sqrt(k^2+n/v1);
   q2 = sqrt(k^2+n/v2);
   T = 0;
   g = 9.81;
   
   r = -((g*k/n^2)*((alpha1-alpha2)+k^2*T/(g*(rho1+rho2)))+1)*(alpha2*q1 + alpha1*q2 - k) - 4*k*alpha1*alpha2 + (4*k^2/n)*(alpha1*v1 - alpha2*v2)*(alpha2*q1-alpha1*q2+k*(alpha1-alpha2)) + (4*k^3/n^2)*(alpha1*v1-alpha2*v2)^2*(q1-k)*(q2-k);
   #r = (k/n^2)*(rho2-rho1)*g-T*k^3/n^2-(2*k^2/n)*(mu2+mu1)-rho2-rho1;
endfunction

n0 = 1.0

#[n,obj,info,iter,nf,lambda] = sqp(n0,@phi)
[n, fval, info] = fsolve (@phi, n0)

   porosity = 0.8;
   k = sqrt(2*4*pi^2);
   rho1 = 1000;
   rho2 = 1000*porosity+8000*(1-porosity);
   mu1 = 0.218375*rho1
   mu2 = mu1*(1-(1-porosity)/0.63)^(-2.5*0.63)
   alpha1 = rho1/(rho1+rho2);
   alpha2 = rho2/(rho1+rho2);
   v1 = mu1/rho1;
   v2 = mu2/rho2;
   q1 = sqrt(k^2+n/v1);
   q2 = sqrt(k^2+n/v2);
   T = 0;
   g = 9.81;
   

   r = -((g*k/n^2)*((alpha1-alpha2)+k^2*T/(g*(rho1+rho2)))+1)*(alpha2*q1 + alpha1*q2 - k) - 4*k*alpha1*alpha2 + (4*k^2/n)*(alpha1*v1 - alpha2*v2)*(alpha2*q1-alpha1*q2+k*(alpha1-alpha2)) + (4*k^3/n^2)*(alpha1*v1-alpha2*v2)^2*(q1-k)*(q2-k)
   #r = (k/n^2)*(rho2-rho1)*g-T*k^3/n^2-(2*k^2/n)*(mu2+mu1)-rho2-rho1

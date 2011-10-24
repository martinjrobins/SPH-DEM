//setup arrays
const double DX = (RMAX[0]-RMIN[0])/NX;
const double DY = (RMAX[1]-RMIN[1])/NY;

Array<NDIM,Cparticle *> pp(NX+1,NY+1);
Array<NDIM,double> s(NX+1,NY+1);
Array<NDIM,double> phi(NX+1,NY+1);

int pn = 0;
for (int i=-2;i<=NX+2;i++) {
   for (int j=-2;j<=NY+2;j++) {
      if ((j>=0)||(j<=NY)||(i>=0)||(i<=NX)) {
         pp(i,j) = &ps[pn];
      }
      pn++;
}

for (int i=0;i<=NX;i++) {
   for (int j=0;j<=NY;j++) {
      if ((j==0)||(j==NY)||(i==0)||(i==NX)) {
         s(i,j) = 0;
      } else {
         s(i,j) = (pp(i-1,j)->v[0]-pp(i+1,j)->v[0])/(2*DX) + (pp(i,j-1)->v[1]-pp(i,j+1)->v[1])/(2*DY);
      }
   }
}

cout << "solving elliptical equation..."<<endl;
for (int iteration;iteration<10;iteration++){
   double maxDelPhi = 0;
   cout << "iteration "<<iteration<<endl;
   for (int i=0;i<=NX;i++) {
      for (int j=0;j<=NY;j++) {
         double oldPhi = phi(i,j);
         if ((i==0)&&(j==0)) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j+1);
         } else if ((i==0)&&(j==NY)) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j-1);
         } else if ((i==NX)&&(j==0)) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j+1);
         } else if ((i==NX)&&(j==NY)) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j-1);
         } else if (i==0) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(2*pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
         } else if (i==NX) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(2*pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
         } else if (j==0) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(pow(DX,2)+2*pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+(pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*phi(i,j+1);
         } else if (j==NY) {
            phi(i,j) = -s(i,j)+(pow(DY,2)/(pow(DX,2)+2*pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+(pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*phi(i,j-1);
         } else {
            phi(i,j) = -0.5*s(i,j)+0.5*(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+0.5*(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
         }
         double delPhi = phi(i,j)-oldPhi;
         if (delPhi > maxDelPhi) maxDelPhi = delPhi;
      }
   }
   cout << "\tmaxDelPhi = "<<maxDelPhi<<endl;
}

for (int i=1;i<NX;i++) {
   for (int j=1;j<NY;j++) {
      vect u;
      u[0] = (phi(i+1,j)-phi(i-1,j))/(2*DX);
      u[1] = (phi(i,j+1)-phi(i,j-1))/(2*DY);
      pp->v -= u;
   }
}

for (int i=0;i<=NX;i++) {
   for (int j=0;j<=NY;j++) {
      if ((j==0)||(j==NY)||(i==0)||(i==NX)) {
         s(i,j) = 0;
      } else {
         s(i,j) = (pp(i-1,j)->v[0]-pp(i+1,j)->v[0])/(2*DX) + (pp(i,j-1)->v[1]-pp(i,j+1)->v[1])/(2*DY);
      }
   }
}

cout << "maximum -div v = "<<max(abs(s))<<endl;      

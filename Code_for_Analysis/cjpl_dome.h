#include "math.h"

#ifndef INCLUDE_CONSTANT_H
#define INCLUDE_CONSTANT_H
#include "constant.h"
#endif

//const double Lx = 18.0*100000.0;// cm
//const double Ly = 100.0*100000.0;// cm
//const double H = 2.4*100000.0;// cm
const double Lx = 18.0; //km
const double Ly = 100.0; //km
const double H = 2.4; //km
//
double Rx = ((H*H)+(Lx*Lx/4.0))/(2.0*H);
double Ry = ((H*H)+(Ly*Ly/4.0))/(2.0*H);
//

double domehigh(double x, double y)
{
double z;

if((abs(x)<Lx/2.0)&&(abs(y)<Ly/2.0)&&(sqrt(Rx*Rx-x*x)>(Rx-H))&&(sqrt(Ry*Ry-y*y)>(Ry-H)))
{ z = ( sqrt(Rx*Rx-x*x)-(Rx-H) )*( sqrt(Ry*Ry-y*y)-(Ry-H) )/H; }
else { z = 0.0; }

if(z<0.0) { z = 0.0; }

//if((x>(Ly/2.0-1.0))&&(y<(-Ly/2.0+1.0))) { z = 15; }

return z;
}

double Lcjpl(double x_dir, double y_dir, double z_dir)
{
// x/x_dir = y/y_dir = z/z_dir

double x_solution, y_solution, z_solution;
double L_cjpl;

if(z_dir<0.0) { z_dir = abs(z_dir); }
if(z_dir<1e-16) { z_dir = 1e-16; }

if(x_dir<0.0) { x_dir = abs(x_dir); }
if(y_dir<0.0) { y_dir = abs(y_dir); }

x_solution = 0.0;
y_solution = 0.0;
z_solution = H;
L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);

if((abs(z_dir)>=abs(x_dir))&&(abs(z_dir)>=abs(y_dir)))
{
  TF1 *fz = new TF1("fz",Form("domehigh(%e*x, %e*x) - x", (x_dir/z_dir), (y_dir/z_dir)),0,H);
  z_solution = fz->GetX(0.0);
  x_solution = z_solution*(x_dir/z_dir);
  y_solution = z_solution*(y_dir/z_dir);
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  fz->Delete();
  //printf("Zsol ");
}
else if((abs(x_dir)>=abs(y_dir))&&(abs(x_dir)>=abs(z_dir)))
{
  TF1 *fx = new TF1("fx",Form("domehigh(x, %e*x) - %e*x", (y_dir/x_dir), (z_dir/x_dir)),0,Lx);
  x_solution = fx->GetX(0.0);
  y_solution = x_solution*(y_dir/x_dir);
  z_solution = x_solution*(z_dir/x_dir); 
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  fx->Delete();
  //printf("Xsol "); 
}
else if((abs(y_dir)>=abs(x_dir))&&(abs(y_dir)>=abs(z_dir)))
{
  TF1 *fy = new TF1("fy",Form("domehigh(%e*x, x) - %e*x", (x_dir/y_dir), (z_dir/y_dir)),0,Ly);
  y_solution = fy->GetX(0.0);
  x_solution = y_solution*(x_dir/y_dir);
  z_solution = y_solution*(z_dir/y_dir);
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  fy->Delete();
  //printf("Ysol ");
}
else 
{
  x_solution = 0.0;
  y_solution = 0.0;
  z_solution = H;
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  //printf("Nsol ");
}

//printf("%e %e %e L=%e\n",x_solution,y_solution,z_solution,L_cjpl);
//
return 100000.0*L_cjpl;
}

double Lcjpl_NonNL(double x_dir, double y_dir, double z_dir)
{
// x/x_dir = y/y_dir = z/z_dir
double NFFCJPL = sqrt(x_dir*x_dir+y_dir*y_dir+z_dir*z_dir);
x_dir = x_dir/NFFCJPL;y_dir = y_dir/NFFCJPL;z_dir = z_dir/NFFCJPL;
    
double x_solution, y_solution, z_solution;
double L_cjpl;

if(z_dir<0.0) { z_dir = abs(z_dir); }
if(z_dir<1e-16) { z_dir = 1e-16; }

if(x_dir<0.0) { x_dir = abs(x_dir); }
if(y_dir<0.0) { y_dir = abs(y_dir); }

x_solution = 0.0;
y_solution = 0.0;
z_solution = H;
L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);

if((abs(z_dir)>=abs(x_dir))&&(abs(z_dir)>=abs(y_dir)))
{
  TF1 *fz = new TF1("fz",Form("domehigh(%e*x, %e*x) - x", (x_dir/z_dir), (y_dir/z_dir)),0,H);
  z_solution = fz->GetX(0.0);
  x_solution = z_solution*(x_dir/z_dir);
  y_solution = z_solution*(y_dir/z_dir);
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  fz->Delete();
  //printf("Zsol ");
}
else if((abs(x_dir)>=abs(y_dir))&&(abs(x_dir)>=abs(z_dir)))
{
  TF1 *fx = new TF1("fx",Form("domehigh(x, %e*x) - %e*x", (y_dir/x_dir), (z_dir/x_dir)),0,Lx);
  x_solution = fx->GetX(0.0);
  y_solution = x_solution*(y_dir/x_dir);
  z_solution = x_solution*(z_dir/x_dir);
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  fx->Delete();
  //printf("Xsol ");
}
else if((abs(y_dir)>=abs(x_dir))&&(abs(y_dir)>=abs(z_dir)))
{
  TF1 *fy = new TF1("fy",Form("domehigh(%e*x, x) - %e*x", (x_dir/y_dir), (z_dir/y_dir)),0,Ly);
  y_solution = fy->GetX(0.0);
  x_solution = y_solution*(x_dir/y_dir);
  z_solution = y_solution*(z_dir/y_dir);
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  fy->Delete();
  //printf("Ysol ");
}
else
{
  x_solution = 0.0;
  y_solution = 0.0;
  z_solution = H;
  L_cjpl = sqrt(x_solution*x_solution+y_solution*y_solution+z_solution*z_solution);
  //printf("Nsol ");
}

    return L_cjpl;//km
}


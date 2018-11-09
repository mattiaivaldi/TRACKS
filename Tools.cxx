#include "Tools.h"
#include "TMath.h"

using namespace TMath;

double *hit_point(double x0, double y0, double z0, double theta, double phi, double R) {
    
    static double hit[3];
    double c1 = Sin(theta)*Cos(phi), c2 = Sin(theta)*Sin(phi), c3 = Cos(theta); //direction cosines 
    double delta = 2*x0*y0*c1*c2 - c1*c1*y0*y0 + c1*c1*R*R -c2*c2*x0*x0 + c2*c2*R*R; //delta of II degree equation ( >= 0 by construction)
    double t_p = (-(x0*c1 - y0*c2) + Sqrt(delta))/(c1*c1 + c2*c2); //solution with the "+" sign 
    double t_m = (-(x0*c1 - y0*c2) - Sqrt(delta))/(c1*c1 + c2*c2); //solution with the "-" sign 

    //calculate the values of the intersection points (x,y,z)
    if (t_p >= 0) {
        hit[0] = x0 + c1*t_p;
        hit[1] = y0 + c2*t_p;
        hit[2] = z0 + c3*t_p;
    }
    
    else {
        hit[0] = x0 + c1*t_m;
        hit[1] = y0 + c2*t_m;
        hit[2] = z0 + c3*t_m;
    }

    //returns the pointer to the first element of the array hit which contains the values x, y and z
    return hit;
}

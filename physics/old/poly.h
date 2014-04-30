#include <string>


// Evaluate polynomials
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// using Horner's rule, ascending powers
double polyval(double c[], int n, const double x){
   p = c[n]; for(int i=n-1; i>=0; i--){p = c[i] + x*p;}
}
// useful overloads
void polyval(double c[], int n, const double x, double &p){p = polyval(c,n,x);}
double polyval(double c[], int n, const double x, String direction){
   if strcmp(direction,"ascending")
      return polyval(c,n,x);
   else if strcmp(direction,"descending"){
      p = c[0]; for(int i=1; i<n; i++){p = c[i] + x*p;}
      return p;
   }
}
void polyval(double c[], int n, const double x, double &p, String direction)
     {p = polyval(c,n,x,direction);}



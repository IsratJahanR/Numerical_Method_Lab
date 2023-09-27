///bisection
/*
#include<bits/stdc++.h>
using namespace std;

double fun(double x)
{
   return x*x*x+2*x*x-x+8;
}
int main()
{

   double xl=-100,xu=100,xr;
   double fxl,fxu,fxr;

   int it=1;

   fxl = fun(xl);
   fxu = fun(xu);

   while(it<=100)
   {
        xr =(xl+xu)/2.0;
      cout<<" "<<it<<"       "<<xr<<endl;

      fxr = fun(xr);
      if(fxl*fxr<0)
      {
          xu=xr;
          fxu=fxr;
      }
      else
      {
          xl=xr;
          fxl=fxr;
      }

      it++;

   }

}

*/
///false position:
/*
double getfun(double x)
{
  double ans = (667.38/x)*(1-exp(-0.146843*x))-40;
  return ans;
}

int main()
{

   double xl=12,xu=16,xr=0,es=0.5,ea=100,et=0.0;

   double fxl,fxu,fxr;

   fxl=getfun(xl);
   fxu=getfun(xu);


   int it=1;

   while(it<=100)
   {

      double  xr1 =xu-((fxu*(xl-xu))/(fxl-fxu));

      xr=xr1;
      cout<<it<<"     "<<xr<<endl;




      fxr=getfun(xr);
      if(fxl*fxr<0)
      {
          xu=xr;
          fxu=fxr;
      }
      else
      {
          xl=xr;
          fxl=fxr;
      }

      it++;


   }

}
*/
///Sacand Method:
/*
#include<bits/stdc++.h>
using namespace std;

double getfun(double x)
{
  double ans = exp(-x)-x;
  return ans;
}

int main()
{

   double xi1=0,xi2,es=1.94,ea=100,et=0.0,xt;

   xi1=0;
   xi2=1.0;


   vector<double>x;
   x.push_back(xi1);
   x.push_back(xi2);
   int it=1;
   while(it<5)
   {

      double  xi =x[it]-(getfun(x[it])*(x[it-1]-x[it]))/(getfun(x[it-1])-getfun(x[it]));
      x.push_back(xi);
      cout<<" "<<it<<"       "<<x[it+1]<<endl;
      it++;

   }

}
*/

///trapizoidal
/*
double getfun(double x)
{
  double ans = 0.2+25*x-200*x*x+675*x*x*x-900*x*x*x*x+400*x*x*x*x*x;

  return ans;
}

int main()
{

   double a=0,b=0.8;
   double sum=getfun(a)+getfun(b);
   double n=2;
   double h= (b-a)/n;

   double strt=h;
   while(strt<b)
   {
       sum+=2*getfun(strt);
       strt+=h;
   }


   double ans =h*(sum/2);
   cout<<ans<<endl;


}*/

///simpson
/*
double getfun(double x)
{
  double ans = 0.2+25*x-200*x*x+675*x*x*x-900*x*x*x*x+400*x*x*x*x*x;

  return ans;
}
void simpson_1_3()
{
   double a=0,b=0.8;
   double sum=getfun(a)+getfun(b);
   double n=4;
   double h= (b-a)/n;

   double strt=h;
   int it=1;
   while(strt<b)
   {
       if(it%2)sum+=4*getfun(strt);
       else sum+=2*getfun(strt);
       strt+=h;
       it++;
   }


   double ans =h*(sum/3);
   cout<<ans<<endl;
}
void simpson_3_8()
{
   double a=0,b=0.8;
   double sum=getfun(a)+getfun(b);
   double n=3;
   double h= (b-a)/n;

   double strt=h;
   int it=1;
   while(strt<b)
   {
       if(it%3)sum+=3*getfun(strt);
       else sum+=2*getfun(strt);
       strt+=h;
       it++;
   }


   double ans =((3*h)/8)*sum;
   cout<<ans<<endl;
}
void fun_1_3_mul()
{
    double a = 0, b = 2, n = 4; // 4.75
    cin >> a >> b >> n;
    double h = (b-a)/n;
    vector <double> v(100);
    int j = 0;
    for (double i = a; i <= b; i += h)
    {
        v[j++] = getfun(i);
    }
    j = 0;
    double ans = 0;
    for (int i = 0; i < n; i++)
    {
        //ans += h*((v[j]+v[j+1])/2);
        ans += ((v[j]+v[j+1]*4+v[j+2]));
        j += 2;
    }
    cout << ans*((2*h)/6) << endl;
}
int main()
{

simpson_3_8();


}*/


///with error:

/*
#include<bits/stdc++.h>
using namespace std;
void inputFunction()
{
    //cout<<"667.38/x(1-e^(-0.146843*x))-40"<<endl;
}
double getfun(double x)
{
  double ans = exp(-x)-x;
 // if(x==1.0)cout<<ans<<endl;
  return ans;
}

double getEA(double pre,double curr)
{
    double ans =abs((1-pre/curr))*100.0;
    //cout<<ans<<endl;
    return ans;
}
double getET(double xi,double xt)
{
   return abs(1-xi/xt)*100.0;
}
int main()
{
   inputFunction();
   double xi1=0,xi2,es=1.94,ea=100,et=0.0,xt;

   xt=0.567;
   ea = 100.0;

   xi1=0;
   xi2=1.0;
   cout<<"Iter.     Xi         ea%        et% "<<endl;
   cout<<"-------------------------------------------------------"<<endl;

   vector<double>x;
   x.push_back(xi1);
   x.push_back(xi2);
   int it=1;
   while(ea>et)
   {

      double  xi =x[it]-(getfun(x[it])*(x[it-1]-x[it]))/(getfun(x[it-1])-getfun(x[it]));
      //cout<<xi<<endl;
      x.push_back(xi);

      //it++;
      ea=getEA(x[it],x[it+1]);
      et=getET(x[it+1],xt);
      cout<<" "<<it<<"       "<<x[it+1]<<"        "<<" "<<"        "<<et<<endl;

     it++;

   }

}

newton rapspn:
#include<bits/stdc++.h>
using namespace std;
void inputFunction()
{
    //cout<<"667.38/x(1-e^(-0.146843*x))-40"<<endl;
}
double getfun(double x)
{
  double ans = exp(-x)-x;
  return ans;
}
double get_dif_fun(double x)
{
  double ans = -exp(-x)-1;
  return ans;
}

double getEA(double pre,double curr)
{
    double ans =abs((1-pre/curr))*100.0;
    //cout<<ans<<endl;
    return ans;
}
double getET(double xi,double xt)
{
   return abs(1-xi/xt)*100.0;
}
int main()
{
   inputFunction();
   double xi=0,es=1.94,ea=100,et=0.0,xt;

   xt=0.567;
   ea = 100.0;

   int it=1;
   cout<<"Iter.     Xi         ea%        et% "<<endl;
   cout<<"-------------------------------------------------------"<<endl;


   while(ea>es)
   {

      double  xi1 =xi-getfun(xi)/get_dif_fun(xi);


      ea=getEA(xi,xi1);
      xi=xi1;
      et=getET(xi,xt);
      cout<<" "<<it<<"       "<<xi<<"        "<<ea<<"        "<<et<<endl;

      it++;

   }

}


//modified sacand



#include<bits/stdc++.h>
using namespace std;
void inputFunction()
{
    //cout<<"667.38/x(1-e^(-0.146843*x))-40"<<endl;
}
double getfun(double x)
{
  double ans = exp(-x)-x;
  return ans;
}

double getEA(double pre,double curr)
{
    double ans =abs((1-pre/curr))*100.0;
    return ans;
}
double getET(double xi,double xt)
{
   return abs(1-xi/xt)*100.0;
}
int main()
{
   inputFunction();
   double xi=0,es=1.94,ea=100,et=0.0,xt;
   double del=0.01;
   xt=0.567;
   ea = 100.0;

   cout<<"Iter.     Xi         ea%        et% "<<endl;
   cout<<"-------------------------------------------------------"<<endl;

   vector<double>x;
   x.push_back(xi);
   int it=0;
  // cout<<x[it]<<endl;
   while(ea>es)
   {

      double  xi =x[it]-(del*x[it]*getfun(x[it]))/(getfun(x[it]+del*x[it])-getfun(x[it]));
      cout<<getfun(x[it]+del*x[it])-getfun(x[it])<<endl;
      x.push_back(xi);

      it++;
      ea=getEA(x[it-1],x[it]);
      et=getET(x[it],xt);
      cout<<" "<<it<<"       "<<x[it]<<"        "<<" "<<"        "<<et<<endl;

   }

}


///fixed point iteration

#include<bits/stdc++.h>
using namespace std;
void inputFunction()
{
    //cout<<"667.38/x(1-e^(-0.146843*x))-40"<<endl;
}
double getfun(double x)
{
  double ans = exp(-x);
  return ans;
}
double getEA(double pre,double curr)
{
    double ans =abs((1-pre/curr))*100.0;
    //cout<<ans<<endl;
    return ans;
}
double getET(double xi,double xt)
{
   return abs(1-xi/xt)*100.0;
}
int main()
{
   inputFunction();
   double xi=0,es=1.94,ea=100,et=0.0,xt;
   //cin>>xi;
   //cin>>es;
   xt=0.567;
   ea = 100.0;

   int it=2;
   cout<<"Iter.     Xi         ea%        et% "<<endl;
   cout<<"-------------------------------------------------------"<<endl;


   while(ea>es)
   {

      double  xi1 =getfun(xi);


      ea=getEA(xi,xi1);
      xi=xi1;
      et=getET(xi,xt);
      cout<<" "<<it<<"       "<<xi<<"        "<<ea<<"        "<<et<<endl;

      it++;

   }

}

//false position

#include<bits/stdc++.h>
using namespace std;
void inputFunction()
{
    cout<<"667.38/x(1-e^(-0.146843*x))-40"<<endl;
}
double getfun(double x)
{
  double ans = (667.38/x)*(1-exp(-0.146843*x))-40;
  return ans;
}
double getEA(double pre,double curr)
{
    double ans =abs((1-pre/curr))*100.0;
    //cout<<ans<<endl;
    return ans;
}
double getET()
{

}
int main()
{
   inputFunction();
   double xl=12,xu=16,xr=0,es=0.5,ea=100,et=0.0;
   //cin>>xl>>xu;
   //cin>>es;
   ea = 100.0;


   double fxl,fxu,fxr;

   fxl=getfun(xl);
   fxu=getfun(xu);

   //xr =xu-((fxu*(xl-xu))/(fxl-fxu));
   int it=1;
   cout<<"Iter.     Xl       Xu        Xr         ea%        et% "<<endl;
   cout<<"-------------------------------------------------------"<<endl;

   //cout<<exp(1)<<endl;
   while(ea>es)
   {

      double  xr1 =xu-((fxu*(xl-xu))/(fxl-fxu));
      ea = getEA(xr,xr1);
      xr=xr1;
      cout<<" "<<it<<"       "<<xl<<"        "<<xu<<"        "<<xr<<"        "<<ea<<"        "<<et<<endl;




      fxr=getfun(xr);
      if(fxl*fxr<0)
      {
          xu=xr;
          fxu=fxr;
      }
      else
      {
          xl=xr;
          fxl=fxr;
      }

      it++;


   }

}



///

#include<bits/stdc++.h>
using namespace std;
void inputFunction()
{
    cout<<"667.38/x(1-e^(-0.146843*x))-40"<<endl;
}
double getfun(double x)
{
  double ans = (667.38/x)*(1-exp(-0.146843*x))-40;
  return ans;
}
double getEA(double xl,double xu)
{
    double ans =((xu-xl)/(xu+xl))*100.0;
    //cout<<ans<<endl;
    return ans;
}
double getET()
{

}
int main()
{
   inputFunction();
   double xl=12,xu=16,xr,es=0.5,ea=100,et=0.0;
   //cin>>xl>>xu;
   //cin>>es;
   ea = getEA(xl,xu);


   double fxl,fxu,fxr;

   fxl=getfun(xl);
   fxu=getfun(xu);


   int it=1;
   cout<<"Iter.     Xl       Xu        Xr         ea%        et% "<<endl;
   cout<<"-------------------------------------------------------"<<endl;

   //cout<<exp(1)<<endl;
   while(ea>es)
   {
        xr =(xl+xu)/2.0;
        ea = getEA(xl,xu);
      cout<<" "<<it<<"       "<<xl<<"        "<<xu<<"        "<<xr<<"        "<<ea<<"        "<<et<<endl;

      fxr=getfun(xr);
      if(fxl*fxr<0)
      {
          xu=xr;
          fxu=fxr;
      }
      else
      {
          xl=xr;
          fxl=fxr;
      }

      it++;

   }

}
*/

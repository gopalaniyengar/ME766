#include<omp.h>
#include<math.h>
#include<iostream>
#include<chrono>

#define a -3.14159265359/2
#define b +3.14159265359/2

using namespace std;


double f(double x)
{
	double k = cos(x);
	return k;
}

double trap(double x1, double x2)
{
	double val=0;
	val=((f(x1)+f(x2))*(b-a))/2;
	return val;
}

double monte(double x1, double x2, double x3, double x4)
{
	double val=0;
	val=((f(x1)+f(x2)+f(x3)+f(x4))*(b-a))/4;
	return val;
}


double calc_trap1(int n)
{
	double h=(double)(b-a)/n;
	double sum=(f(a)+f(b))/2;
	double xiter,fx;
	#pragma omp parallel for default(none), private(xiter,fx), shared(n,h), reduction(+:sum)
	for(int i=1; i<n;i++)
	{
		xiter=(double)a+(h*i);
		fx=(double)cos(xiter);
		sum+=(double)(h*fx);		
	}
	return sum;
}

/*double calc_trap2(int n)
{
	double h=(double)(b-a)/n;
	double sum=0.0;
	double X[n+1];
	#pragma omp parallel for default(none), shared(X,n,h)
	for(int i=0;i<=n;i++)
	{
		X[i]=a+(i*h);
	}
	#pragma omp parallel for default(none), shared(n,X), reduction(+:sum)
	for(int i=0; i<n; i++)
	{
		sum+=trap(X[i],X[i+1]);
	}
	sum=sum/n;
	return sum;
}*/


double calc_monte(int n)
{
	//if(n<=3)
	//{n=3;}
	//else if(n>3 && (n%3)<=1)
	//{n-=(n%3);}
	//else if(n>3 && (n%3)>1)
	//{n+=1;}			//reasonable approximation of n to multiple of 3(for application of composite Monte-Carlo method)
	
	double sum=0.0;
	double h=(double)(b-a)/n; 

	double xiter,r,fx,random;

	#pragma omp parallel for default(none), private(xiter,fx,random,r), shared(n,h), reduction(+:sum)
	for(int i=0;i<=n;i++)
	{
		float random = ((float) rand()) / (float) RAND_MAX;
   		r = random * (b-a);		
		xiter=a+r;
		fx=cos(xiter);
		sum+=(h*fx);	
	}
	return sum;
}

int main()
{
    double traptimes[7][4];
    double montetimes[7][4];	
    double tenindex[7];
    
    for(int j=3;j<=9;j++)
	{
		int k=j-3;
		long long int ndiv=pow(10,j);
		cout<<"POWER OF 10: "<<j<<endl;
		for(int threads=2;threads<=8;threads+=2)
		{
			omp_set_num_threads(threads);
			double t1=omp_get_wtime()*1000;
			double integral1=calc_trap1(ndiv);
			double t2=omp_get_wtime()*1000;
			double t3=omp_get_wtime()*1000;
			double integral2=calc_monte(ndiv);
			double t4=omp_get_wtime()*1000;
			double timetrap=abs(t2-t1);
			double timemonte=abs(t4-t3);
		cout<<"Trapezoid: "<<integral1<<", time: "<<timetrap<<"ms"<<endl;
		cout<<"Monte Carlo: "<<integral2<<", time: "<<timemonte<<"ms"<<endl;
		int m=(threads/2)-1;
		traptimes[k][m]=timetrap;
		montetimes[k][m]=timemonte;
		}
		cout<<endl<<endl;
	}
	
    omp_set_num_threads(8);
	for(long int j=1; j<=1000000;j*=10)
	{
	cout<<j<<": ";
	double integral1= calc_trap1(j);
	double integral2= calc_monte(j);
	cout<<"Trapezoidal= "<<integral1<<" ;  Monte Carlo= "<<integral2<<endl;
	}
} 
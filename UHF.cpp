#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <iostream>
#include <random>
#include <time.h>

using namespace std;
using namespace NTL;

int main()
{
clock_t tbegin, tend;
tbegin=clock();
long n=4160;
long err=80;
ZZ q1;
GenPrime(q1, n, err);
cout<<"q1="<<q1<<"\n";
//randomly choose a number
ZZ_p::init(q1);
mat_ZZ_p x;
random(x,1,1);
mat_ZZ_p a;
random(a,1,1);
mat_ZZ_p b;
random(b,1,1);

cout<<"x="<<x<<"\n";
//transform into binary and use to encrypt
ZZ temp_x;
temp_x=conv<ZZ>(x[0][0]);
mat_ZZ_p q1_2;
q1_2.SetDims(1,n);
for(int i=0;i<n;i++)
	{
		q1_2[0][i]=temp_x%2;
		if (q1_2[0][i]!=0)
			{
				temp_x=(temp_x-1)/2;
			}
			else
				{
					temp_x=temp_x/2;
				}
	}
cout<<"q1_2="<<q1_2<<"\n";//obtain a string 



long n_y=4160*5/130-2*13;
cout<<"n_y="<<n_y<<"\n";
ZZ q2;
GenPrime(q2, n_y, err);
cout<<"q2="<<q2<<"\n";

ZZ_p::init(q2);
ZZ_p y_hash;
y_hash=conv<ZZ_p>(conv<ZZ>(a[0][0]*x[0][0])+conv<ZZ>(b[0][0]));
cout<<"y_hash="<<y_hash<<"\n";
//transform into binary
ZZ temp_y;
temp_y=conv<ZZ>(y_hash);
mat_ZZ_p q2_2;
q2_2.SetDims(1,n_y);
for(int i=0;i<n_y;i++)
	{
		q2_2[0][i]=temp_y%2;
		if (q2_2[0][i]!=0)
			{
				temp_y=(temp_y-1)/2;
			}
			else
				{
					temp_y=temp_y/2;
				}
	}
cout<<"q2_2="<<q2_2<<"\n";//obtain a string 
tend=clock();
double t_3;
t_3=double(tend-tbegin)/CLOCKS_PER_SEC;
cout<<"t_3="<<t_3<<"\n";

}
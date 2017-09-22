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

//分享陷门函数,用户输入身份，函数输出用户的私钥分量，最好是输出整个私钥，但是这个和参数w相关（可以写两个函数共同完成该功能）
//需要定义系数矩阵，输入变量（身份）,多项式的最高次幂（T-1）
//计算的陷门矩阵是d*w
//下面求解各用户的分享陷门矩阵
//算法是应该用w个d*T的的矩阵分别乘上T*1的矩阵，得到的w个d*1的矩阵组成用户的分享陷门矩阵；
//所以这里要先分割矩阵，然后再乘T*1的矩阵
mat_ZZ_p sharefun1(mat_ZZ_p coff_share, long T_share, long id_share, long d_share, long w_share)
	{
		int i, j, k, kk, kkk;
		mat_ZZ_p id_matrix;//定义身份向量T_share*1；
		id_matrix.SetDims(T_share,1);

		for(i=0; i<T_share; i++)
			{
				id_matrix[i][0]=power(conv<ZZ_p>(conv<ZZ>(id_share)), i);
			}
		cout<<id_matrix<<"\n";

		mat_ZZ_p td_share;
		td_share.SetDims(d_share, w_share);//定义用户分享陷门矩阵d*w
		//分割矩阵，分别相乘
		for(i=0; i<w_share; i++)
			{
				kk=i*T_share;
				kkk=(i+1)*T_share;
				for(j=0;j<d_share;j++)
					{
						for(k=kk;k<kkk;k++)
							{
								td_share[j][i]=td_share[j][i]+coff_share[j][k]*id_matrix[k-i*T_share][0];
							}
					}
			}
		return td_share;
	}


//为了方便，直接先生成w个随机的系数矩阵，因为每个随机矩阵的第一列构成了陷门，而且都是Z_q上的随机选择

int main()
{
clock_t tbegin, tend;
int i,j,k; 
long n, d, w, m, p, N;
long T, eta;//n代表加密用的随机串的长度；d代表安全参数，即私钥矩阵的行数；w代表私钥矩阵的列数；
//m代表m*l=n；p消息矩阵的模数；N用户总数；T解密至少需要的用户数； eta代表N^2(error) (N!)^2；
ZZ q;

n=1100, d=256, p=2027;
w=m=100;

q=1324909189, N=3, T=3;
int eta2=1;
for (i=1;i<=N;i++)
	{
		eta2=eta2*i;
	}
eta=eta2*eta2;

ZZ_p::init(q);
mat_ZZ_p A, Z;

//mat_ZZ_p BB;
//random(BB, d, T*w);
//(pk) random A and (sk) random Z 
tbegin=clock();
random(A, n, d);
cout<<"A="<<A<<"\n";
Z.SetDims(d, w);
mat_ZZ_p BB;//定义BB用来储存d行w*t列的随机矩阵，其中第i列是陷门矩阵i=1,...,w-1;因为很难在ntl中定义三维矩阵，只能将w个矩阵链接起来，
//使用时在分割出来
random(BB, d, T*w);
cout<<"BB="<<BB<<"\n";
for (i=0; i<w; i++){for(j=0;j<d;j++){Z[j][i]=BB[j][i*T];}}
cout<<"Z="<<Z<<"\n";
//思考生成随机矩阵的数量随输入w的变化而变化

//定义N个用户的私钥矩阵
mat_ZZ_p ZZ0, ZZ1, ZZ2, ZZ3, ZZ4;

long id0, id1, id2, id3, id4;
id0=1;id1=2;id2=3;id3=4;id4=5;
//根据share函数求分享的矩阵
//分享陷门，计算N个用户各自的私钥ZZ1,ZZ2,ZZ3,...,ZZN；
ZZ0=sharefun1(BB, T, id0, d, w);
ZZ1=sharefun1(BB, T, id1, d, w);
ZZ2=sharefun1(BB, T, id2, d, w);
ZZ3=sharefun1(BB, T, id3, d, w);
ZZ4=sharefun1(BB, T, id4, d, w);

cout<<"ZZ0="<<ZZ0<<"\n";
cout<<"ZZ1="<<ZZ1<<"\n";
cout<<"ZZ2="<<ZZ2<<"\n";
cout<<"ZZ3="<<ZZ3<<"\n";
cout<<"ZZ4="<<ZZ4<<"\n";
 
mat_ZZ_p A_Z;//def the product of A and Z
mul(A_Z, A, Z);
cout<<"A_Z="<< A_Z <<"\n";

//def the message matrix B
long l=n/m;//n=m*l
//1,def vector b
Mat<RR> b;
b.SetDims(l, 1);
for (i=0;i<l;i++)
	{
		b[i][0]=pow(2, i);
	}
cout<<"b="<<b<<"\n";

mat_RR I;
mat_ZZ_p B;//def message matrix
ident(I,m);
//Mat<RR> B;
B.SetDims(n, m);

int ii=0;
for (j=0; j<m; j++)
	{
		for (i=0; i<m; i++)
			{
				for (k=0;k<l;k++)
					{
						B[ii][j]=conv<ZZ_p>(conv<ZZ>(I[i][j]*b[k][0]*conv<RR>(q)/conv<RR>(p)+0.5));
						ii=ii+1;
						if (ii==n)
							{
								ii=0;
							}
					}
			}
	}
cout<<"B="<<B<<"\n";

//def error matrix
double delta;
delta=1.0/(32*p*n);
unsigned seed=time(NULL);
default_random_engine generator (seed);
normal_distribution<double>distribution(0.0,delta);

Mat<ZZ_p> E;
E.SetDims(n, m);
for (i=0;i<n;i++)
	{
		for (j=0;j<m;j++)
			{
				E[i][j]=conv<ZZ_p>(ZZ(1324909189*distribution(generator)+0.5));
			}
	}

cout <<"E="<<E << "\n";

//def ciphertex matrix C (function index injective)
Mat<ZZ_p> C;
C=A_Z+B+E;
cout<<"C="<<C<<"\n";

tend=clock();
double t_1;
t_1=double(tend-tbegin)/CLOCKS_PER_SEC;//the time is cost by generating the pk sk and shared sk 
//*************************************************
//eva y=x*C
//first randomly choose 0/1 strings
tbegin=clock();
ZZ binary;
binary=2;
ZZ_p::init(binary);
mat_ZZ_p x;
random(x, 1, n);
cout<<"x="<<x<<"\n";

ZZ_p::init(q);
mat_ZZ_p y1, y2;
y2=x*C;
cout<<"y2=x*C"<<y2<<"\n";
y1=x*A;
cout<<"y1=x*A"<<y1<<"\n";
tend=clock();
double t_2;
t_2=double(tend-tbegin)/CLOCKS_PER_SEC;
//*************************************************
//*************************************************
//*************************************************
//def inversion algorithm
tbegin=clock();
ZZ_p eta1;
eta1=inv(conv<ZZ_p>(conv<ZZ>(eta)));
mat_ZZ_p aa;
aa=eta1*y1;
cout<<"aa=eta1*y1"<<aa<<"\n";
//下面计算各个用户的求逆分享
mat_ZZ_p inv0, inv1, inv2, inv3, inv4;
mat_ZZ_p inv0_transpose, inv1_transpose, inv2_transpose, inv3_transpose, inv4_transpose;
mat_ZZ_p ee0, ee1,ee2,ee3,ee4;

ee0.SetDims(1, w);
for(i=0;i<w;i++)
	{
		ee0[0][i]=conv<ZZ_p>(ZZ(1324909189*distribution(generator)+0.5));
	}
cout <<"ee0="<<ee0<<"\n";

inv0_transpose=aa*ZZ0+ee0;//行向量
inv0=transpose(inv0_transpose);
cout <<"inv0_transpose=aa*ZZ0+ee0="<<inv0_transpose<<"\n";//列向量
//inv0_transpose=aa*ZZ0;
//cout <<"inv0_transpose=aa*ZZ0"<<inv0_transpose<<"\n";//列向量
tend=clock();
double t_3;
t_3=double(tend-tbegin)/CLOCKS_PER_SEC;

ee1.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee1[0][i]=conv<ZZ_p>(ZZ(1324909189*distribution(generator)+0.5));
	}
cout <<"ee1="<<ee1<<"\n";

inv1_transpose=aa*ZZ1+ee1;//行向量
inv1=transpose(inv1_transpose);
cout <<"inv1="<<inv1<<"\n";//列向量

ee2.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee2[0][i]=conv<ZZ_p>(ZZ(1324909189*distribution(generator)+0.5));
	}
cout <<"ee2="<<ee2<<"\n";

inv2_transpose=aa*ZZ2+ee2;//行向量
inv2=transpose(inv2_transpose);
cout <<"inv2="<<inv2<<"\n";//列向量

ee3.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee3[0][i]=conv<ZZ_p>(ZZ(1324909189*distribution(generator)+0.5));
	}
cout <<"ee3="<<ee3<<"\n";

inv3_transpose=aa*ZZ3+ee3;//行向量
inv3=transpose(inv3_transpose);
cout <<"inv3="<<inv3<<"\n";//列向量

ee4.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee4[0][i]=conv<ZZ_p>(ZZ(1324909189*distribution(generator)+0.5));
	}
cout <<"ee4="<<ee4<<"\n";

inv4_transpose=aa*ZZ4+ee4;//行向量
inv4=transpose(inv4_transpose);
cout <<"inv4="<<inv4<<"\n";//列向量

//下面设计组合算法combine
//首先输入门限（T）个用户的身份和解密分享（因为T=3,假设输入的是T0=1,T1=2,T2=3；inv0,inv1,inv2）s
//输入用户身份；
tbegin=clock();
long T0=1,T1=2,T2=3;
mat_ZZ_p ID;//定义一个身份矩阵
ID.SetDims(1, T); 
ID[0][0]=T0;
ID[0][1]=T1;
ID[0][2]=T2;

mat_ZZ_p coff_com;//定义一个系数矩阵
coff_com.SetDims(1, T); 

for(i=0;i<T;i++)
	{
		coff_com[0][i]=1;
		for(j=0;j<T;j++)
			{
				if (i==j) 
					{
						continue;
					}
				coff_com[0][i]=coff_com[0][i]*(-ID[0][j])*inv(ID[0][i]-ID[0][j]);
				//此处注意因为ZZ_p类型下ZZ_p a； a=-3; 输出的a=3,并非p-3；一般需要先把a转换成ZZ类型，进行加减，然后？？？？转换成ZZ_p类型因为
				//ZZ a; a=-3; ZZ_p b; b=conv<ZZ_p>(a); 输出b=p-3;？？？？？？似乎没有问题这里，因为在abotdf中最后使用高斯消元法解方程求v时，出现
				//这个问题；
			}
	}
cout<<"coff_com="<<coff_com<<"\n";
cout<<"eta="<<eta<<"\n";
//由于上面用户的身份是1 2 3，下面也用inv0，inv1，inv2；
//将用户的解密分享组合出结果；
mat_ZZ_p yy;//定义新的矩阵用于存放组合后的密文；
yy=eta*coff_com[0][0]*inv0_transpose+eta*coff_com[0][1]*inv1_transpose+eta*coff_com[0][2]*inv2_transpose;
cout<<"yy=(eta*coff_com[0][0])*inv0_transpose"<<yy<<"\n";

mat_ZZ_p Z2;

Z2=y2-yy;//相当于密文向量减去A和Z的内积以及部分噪声；
cout<<"Z2=y2-yy"<<Z2<<"\n";

mat_RR t;
t.SetDims(1, m);//此处应该是w，但是m=w，也就不区分了；

for (i=0;i<m;i++)
	{
		t[0][i]=conv<RR>(conv<ZZ>(Z2[0][i]))/conv<RR>(q);
	}
cout<<"t="<<t<<"\n";

mat_RR closest, xx;
closest.SetDims(1, p);
xx.SetDims(1, m);

//ZZ temp;

for (j=0;j<m;j++)
	{
		for (i=0;i<p;i++)
			{
				closest[0][i]=abs(abs(t[0][j]-double(i)/conv<RR>(p))-0.5);
			}
//cout<<"closest="<<closest<<"\n";
RR max;
max=0;
int mm;
for (i=0; i<p-1;i++)
	{
		if (closest[0][i]>max)
			{
				max=closest[0][i];
				mm=i;
			}
	}
//cout<<"m="<<m<<"\n";
xx[0][j]=to_RR(mm);
}
cout<<"xx="<<xx<<"\n";

//transform into binary
mat_ZZ xxx;
xxx.SetDims(1, n);
ZZ temp;
int iii=0;
for (j=0; j<m; j++)
	{
		temp=conv<ZZ>(xx[0][j]);
		for (i=0;i<l;i++)
			{
				xxx[0][iii]=temp%2;
				if (xxx[0][iii]!=0)
					{
						temp=(temp-1)/2;
					}
					else 
						{
							temp=temp/2;
						}
						iii++;
			}
	}
cout<<"xxx="<<xxx<<"\n";
tend=clock();
double t_4;
t_4=double(tend-tbegin)/CLOCKS_PER_SEC;

//check solution
mat_ZZ solution;
solution.SetDims(1,n);
for (i=0;i<n;i++)
	{
		solution[0][i]=conv<ZZ>(x[0][i])-conv<ZZ>(xxx[0][i]);
	}
long tt=IsZero(solution);
cout<<"solution="<<solution<<"\n";
cout<<"solution is all-zero matrix: "<<tt<<"\n";
cout<<"generate pk and sk: t_1="<<t_1<<"\n";//seconds
cout<<"encrypted time: t_2="<<t_2<<"\n";
cout<<"user's decrypted time: t_3="<<t_3<<"\n";
cout<<"combining time: t_4="<<t_4<<"\n";
//*******************************
//check shamir's scheme
long ttt=IsZero(coff_com[0][0]*ZZ0+coff_com[0][1]*ZZ1+coff_com[0][2]*ZZ2-Z);
cout<<"shamir's scheme is right: "<<ttt<<"\n";

}








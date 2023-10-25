#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
using namespace std;

#define maxBits 8
#define shiftMag 16
#define seed 0
#define ave 100000

static inline uint64_t rdtscp64() {
  uint32_t low, high;
  asm volatile ("rdtscp": "=a" (low), "=d" (high) :: "ecx");
  return (((uint64_t)high) << 32) | low;
}

uint64_t inverse(uint64_t a, uint64_t mod)	//returns the inverse modulus
{
    uint64_t init = mod;
    uint64_t y = 0;
    uint64_t x = 1;
 
    while (a>1)
    {
        uint64_t q = a/mod;
        uint64_t temp = mod;
        mod = a % mod;
        a = temp;
        temp = y;
        y = x - q * y;
        x = temp;
    }
 
    if (x < 0)
        x += init;
 
    return x;
}

vector<uint64_t> precompute(uint64_t n, uint64_t zeta, uint64_t q, int flag)		//used to precompute the roots of unity for zeta
{
	vector<uint64_t> z(n);		
	vector<uint64_t> r{0, n/2};

	for(int i=2; i<n; i*=2)			//bit reverse for n degrees
		for(int j=i-1; j>=0; j--)
			r.push_back(r.at(r.size()-i)+(n/(i*2)));	//add n/(i*2) to the previous (i*2) elements to get the next i*2 elements

	for(int k=0; k<n; k++)
	{
		int temp = 1;

		for(int j=0; j<r[k]; j++)			//z[k] = zeta^brv(k) mod q
			temp = (temp*zeta)%q;	
	
		z[k] = temp;

		if(flag == 1)						//reverse z[k] = zeta^-(brv[k]+1) mod q
			z[k] = inverse((z[k]*zeta)%q, q);
	}

	return z;
}

void print(vector<uint64_t> a)
{
	for(int i=0; i<a.size(); i++)
		cout<<a.at(i)<<" ";
	cout<<endl;
}

vector<uint64_t> generate(uint64_t n, uint64_t q)
{
	vector<uint64_t> t;

	for(int i=0; i<n; i++)
		t.push_back(rand()%q);

	return t;
}

vector<uint64_t> naive(uint64_t n, uint64_t q, vector<uint64_t> z)
{
	vector<uint64_t> p = generate(n,q);
	vector<uint64_t> poly;
	uint64_t temp;
	long long ret = 0; 

	for(int i=0; i<ave; i++)
	{
		poly = p;
		int k = 1;
		int j = 0;
		uint64_t time = rdtscp64();

		for(int l=n/2; l>0; l/=2)                                               //NTT algorithm
			for(int s=0; s<n; s=j+l)
			{
				for(j=s; j<s+l; j++)
				{
					temp = (z[k] * poly[j+l])%q;
					poly[j+l] = (poly[j]+q - temp)%q;
					poly[j] = (poly[j] + temp)%q;
				}

				k++;
			}

		ret += (rdtscp64()-time);
	}

	cout<<"Naive: "<<ret/ave<<endl;
	return poly;
}

vector<uint64_t> unrolled2(uint64_t n, uint64_t q, vector<uint64_t> z)
{
	vector<uint64_t> p = generate(n,q);
	vector<uint64_t> poly;
	uint64_t temp;
	uint64_t ret = 0;

	for(int i=0; i<ave; i++)
	{
		poly = p;
		uint64_t l = n/2;
		uint64_t k = 1;
		uint64_t j = 0;
		uint64_t time = rdtscp64();

		l/=2;
		j = 0;

		for(int s=0; s<n; s=j+l)					//16		32 l =16 s=16+j
		{
			for(j=s; j<s+l; j++)
			{
				temp = (z[k] * poly[j+l])%q;
				poly[j+l] = (poly[j]+q - temp)%q;
				poly[j] = (poly[j] + temp)%q;

			}

			k++;
		}

		l/=2;
		j = 0;

		for(int s=0; s<n; s=j+l)					//8
		{
			for(j=s; j<s+l; j++)
			{
				temp = (z[k] * poly[j+l])%q;
				poly[j+l] = (poly[j]+q - temp)%q;
				poly[j] = (poly[j] + temp)%q;
			}

			k++;
		}

		l/=2;
		j = 0;

		for(int s=0; s<n; s=j+l)					//4
		{
			for(j=s; j<s+l; j++)
			{
				temp = (z[k] * poly[j+l])%q;
				poly[j+l] = (poly[j]+q - temp)%q;
				poly[j] = (poly[j] + temp)%q;
			}

			k++;
		}

		l/=2;
		j = 0;

		for(int s=0; s<n; s+=8)					//2	
		{
			temp = (z[k] * poly[s+2])%q;
			poly[s+2] = (poly[s]+q - temp)%q;
			poly[s] = (poly[s] + temp)%q;

			temp = (z[k] * poly[s+3])%q;
			poly[s+3] = (poly[s+1]+q - temp)%q;
			poly[s+1] = (poly[s+1] + temp)%q;

			k++;

			temp = (z[k] * poly[s+4])%q;
			poly[s+4] = (poly[s+4]+q - temp)%q;
			poly[s+4] = (poly[s+4] + temp)%q;

			temp = (z[k] * poly[s+5])%q;
			poly[s+5] = (poly[s+5]+q - temp)%q;
			poly[s+5] = (poly[s+5] + temp)%q;

			k++;
		}
		

		ret += (rdtscp64()-time);
	}

	cout<<"Unrolled: "<<ret/ave<<endl;
	return poly;
}

int main()	
{
	uint64_t q,zeta,n;		//enter parameters
	n=32;
	q=11777;
	zeta=5558;
	vector<uint64_t> invz = precompute(n, zeta, q, 1);
	vector<uint64_t> z = precompute(n, zeta, q, 0);	

	uint64_t invMod = inverse(n,q);
	uint64_t invq = (1<<shiftMag)/q;

	vector<uint64_t> poly = generate(n,q);	
	
	cout<<"n: "<<n<<" | "<<"Average: "<<ave<<endl;
	naive(n,q,z);
	unrolled2(n,q,z);

	return 0;
}		


//g++ unroll* -O3 -fno-unroll-loops
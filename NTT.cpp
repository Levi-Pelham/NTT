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

vector<uint64_t> forward(uint64_t n, uint64_t q, uint64_t qinv, vector<uint64_t> z)
{
	int k = 1;
	int j = 0; 
	int s = 0;
	vector<uint64_t> poly = {2,5,3,4};
	uint64_t temp;
	uint64_t ret = 0;

	for(int i=0; i<2500; i++)
	{
		uint64_t time = rdtscp64();

		for(int l=n/2; l>0; l/=2)	//NTT algorithm 
		{
			for(int s=0; s<n; s=s+(2*l))	
			{
				for(j=0; j<l; j+=2)
				{
					temp = (z[k] * poly[j+s+l])%q;
					poly[j+s+l] = (poly[j+s]+q - temp);
					poly[j+s] = (poly[j+s] + temp);
				
					temp = (z[k] * poly[j+1+s+l])%q;
					poly[j+1+s+l] = (poly[j+1+s]+q - temp);
					poly[j+1+s] = (poly[j+1+s] + temp);

				}

				k++;
			}
		}

		for(int i=0; i<n; i++)
			poly[i] %= q;

		ret += (rdtscp64()-time);
	}

	cout<<ret/2500<<endl;
	return poly;
}

int main()	
{
	uint64_t q,zeta,n;		//enter parameters
	
	n=2;
	q=17;
	zeta=2;
	vector<uint64_t> invz = precompute(n, zeta, q, 1);
	vector<uint64_t> z = precompute(n, zeta, q, 0);	

	uint64_t invMod = inverse(n,q);
	uint64_t invq = (1<<shiftMag)/q;
	long long total = 0;

	int k, j;
	vector<uint64_t> poly = generate(n,q);
	uint64_t temp;
	uint64_t time;	
	int loop = 0;

	for(int i=0; i<10000; i++)
	{
		int k = 1;
		int j = 0; 
		time = rdtscp64();

		for(int l=n/2; l>0; l/=2)						//NTT algorithm
		{
			for(int s=0; s<n; s+=(2*l))
			{
				for(j=0; j<l; j+=2)
				{
					temp = (z[k] * poly[j+s+l])%q;
					poly[j+s+l] = (poly[j+s]+q - temp)%q;
					poly[j+s] = (poly[j+s] + temp)%q;

					temp = (z[k] * poly[j+s+l+1])%q;
					poly[j+s+l+1] = (poly[j+s+1]+q - temp)%q;
					poly[j+s+1] = (poly[j+s+1] + temp)%q;

					loop++;
					cout<<loop<<" - "<<l<<" : "<<s<<" : "<<j<<endl;
				}

				k++;
			}
		}

		total += rdtscp64()-time;
	}

	cout<<"L = "<<loop<<endl;
	cout<<total/10000<<endl;
	return 0;
}		


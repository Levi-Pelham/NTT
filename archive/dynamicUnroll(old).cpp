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

vector<uint64_t> generate(uint64_t n, uint64_t q)
{
	vector<uint64_t> t;

	for(int i=0; i<n; i++)
		t.push_back(rand()%q);

	return t;
}

void print(vector<uint64_t> a)
{
	for(int i=0; i<a.size(); i++)
		cout<<a.at(i)<<" ";
	cout<<endl;
}

string format(vector<uint64_t> a)
{
	string k = "{";
	for(int i=0; i<a.size(); i++)
	{
		if(i != a.size()-1)
			k += to_string(a.at(i))+",";
		else
			k += to_string(a.at(i))+"}";
	}

	return k;
}

int main()	
{
	uint64_t q,zeta,n,t,index;		//enter parameters
	srand(0);
	// cout<<"n: ";
	// cin>>n;
	// cout<<"q: ";
	// cin>>q;
	// cout<<"zeta: ";
	// cin>>zeta;

	n = 8;
	q = 17;
	zeta = 2;

	vector<uint64_t> z = precompute(n,zeta,q,0);			//precompute all values of zeta and inverse zeta
	vector<uint64_t> invz = precompute(n,zeta,q,1);
	vector<uint64_t> c = generate(n,q);
	uint64_t inv = inverse(n,q);
	uint64_t qinv = (((uint64_t)1<<shiftMag)/q);

	ofstream outfile;			//generate and populate file
	stringstream filename;
	filename<<"NTTu_"<<n<<".cpp";
	outfile.open(filename.str());

	uint64_t l = n/2;	//l loop
	uint64_t factor;	//unrolliing factor
	uint64_t k = 1;

	outfile << "#include <iostream>\nusing namespace std;\n\n static inline uint64_t rdtscp64() {\n\tuint32_t low, high;\n\tasm volatile (\"rdtscp\": \"=a\" (low), \"=d\" (high) :: \"ecx\");\n\treturn (((uint64_t)high) << 32) | low;\n}\n\n";
	outfile << "int main(int argc, char* argv[])\n{\n\tuint64_t z[] = "<<format(z)<<";\n\t";
	outfile << "uint64_t total = 0;";

	outfile << "\n\tfor(int p=0; p<10000; p++)\n\t{\n\n\tuint64_t c[] = "<<format(c)<<";\n\tuint64_t time = rdtscp64();";
	
	outfile << "uint64_t temp;\n\tuint64_t k=1;\n\tuint64_t j=0;\n\n";
	
	while(l>0)		//remove while for single loop testiing <-----
	{
		cout<<l<<" loop factor: ";
		cin>>factor;
		
		if(factor !=1 && factor != n/2)		//factor validation check
			if(factor>l || factor%2 != 0)
			{
				cout<<"invalid factor"<<endl;
				continue;
			}

		uint64_t j=0;		//initialise vairables

		if(factor==1)		//fully rolled
		{
			outfile<<"\tfor(int s=0; s<"<<n<<"; s=j+"<<l<<")\n\t{\n";	

			outfile<<"\t\tfor(j=s; j<s+"<<l<<"; j++)\n";
			outfile<<"\t\t{\n\t\t\ttemp = (z[k] * c[j+"<<l<<"])%"<<q<<";\n\t\t\tc[j+"<<l<<"] = (c[j]+"<<q<<" - temp)%"<<q<<";";
			outfile<<"\n\t\t\tc[j] = (c[j] + temp)%"<<q<<";\n\t\t}\n\t\tk++;\n\t}\n";
		}

		else if(factor == n/2)		//fully unrolled
			for(int s=0; s<n; s=j+l)
			{
				for(j=s; j<s+l; j++)
				{
					outfile<<"\ttemp = (z[k] * c["<<j+l<<"])%"<<q<<";\n";
					outfile<<"\tc["<<j+l<<"] = (c["<<j<<"]+"<<q<<" - temp)%"<<q<<";\n";
					outfile<<"\tc["<<j<<"] = (c["<<j<<"] + temp)%"<<q<<";\n\n";
				}
			
				outfile<<"\tk++;\n";
			}
		
		else if(factor <= l)		//unroll j
		{
			outfile<<"\tfor(int s=0; s<"<<n<<"; s=j+"<<l<<")\n\t{\n";

			//unroll j by factor
			outfile<<"\t\tfor(j=s; j<s+"<<l<<"; j+="<<factor<<")\n\t\t{\n";
			for(j=0; j<factor; j++)
			{
				outfile<<"\t\t\ttemp = (z[k] * c[j+"<<l+j<<"])%"<<q<<";\n";
				outfile<<"\t\t\tc[j+"<<l+j<<"] = (c[j+"<<j<<"] + "<<q<<" - temp)%"<<q<<";\n";
				outfile<<"\t\t\tc[j+"<<j<<"] = (c[j+"<<j<<"] + temp)%"<<q<<";\n";

				j == (factor-1) ? outfile<<"\t\t}\n\t\tk++;\n\t}\n":outfile<<"\n";
			}
		}

		else	//unroll j & s
		{
			//ensure unroll is valid for j
			uint64_t sFact = factor/l;

			outfile<<"\tfor(int s=0; s<"<<n<<"; s=j+"<<l<<")\n\t{\n";

			for(int s=0; s<sFact; s++)
			{
				//update s between j loops
				if(l==1)
					outfile<<"\t\tfor(j=s; j<s+"<<l<<"; j++)\n\t\t{\n";
				
				else
					outfile<<"\t\tfor(j=s; j<s+"<<l<<"; j+="<<l<<")\n\t\t{\n";
				
				for(int j=0; j<l; j++)
				{
					outfile<<"\t\t\ttemp = (z[k] * c[j+"<<l+j<<"])%"<<q<<";\n";
					outfile<<"\t\t\tc[j+"<<l+j<<"] = (c[j+"<<j<<"] + "<<q<<" - temp)%"<<q<<";\n";
					outfile<<"\t\t\tc[j+"<<j<<"] = (c[j+"<<j<<"] + temp)%"<<q<<";\n";

					j == (l-1) ? outfile<<"\t\t}\n\n\t\tk++;\n":outfile<<"\n";
				}

				if(s==sFact-1)
					outfile<<"\t}";
			}
		
		}
		outfile<<"\n";
		l/=2;
	}

	// outfile<<"\tfor(auto x:c)\n\t\tcout<<x<<endl;\n";

	outfile<< "\t\ttotal += rdtscp64()-time;";
	outfile<<"\n\t}\n\tcout<<total/10000<<endl;\n}";

	outfile.close();
	
	return 0;
}		

//g++ unroll* -O3 -fno-unroll-loops
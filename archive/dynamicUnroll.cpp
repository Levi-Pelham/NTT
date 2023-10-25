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
	
	n = 128;
	q = 11777;
	zeta = 5558;

	vector<uint64_t> z = precompute(n,zeta,q,0);			//precompute all values of zeta and inverse zeta
	vector<uint64_t> invz = precompute(n,zeta,q,1);
	vector<uint64_t> c = generate(n,q);
	uint64_t inv = inverse(n,q);
	uint64_t qinv = (((uint64_t)1<<shiftMag)/q);

	ofstream outfile;			//generate and populate file
	stringstream filename;
	filename<<"NTT_"<<n<<".cpp";
	outfile.open(filename.str());

	uint64_t l = n/2;	//l loop
	uint64_t factor;	//unrolliing factor
	uint64_t k = 1;
	uint64_t j=0;

	outfile<<"#include <iostream>\nusing namespace std;\n\n static inline uint64_t rdtscp64() {\n\tuint32_t low, high;\n\tasm volatile (\"rdtscp\": \"=a\" (low), \"=d\" (high) :: \"ecx\");\n\treturn (((uint64_t)high) << 32) | low;\n}\n\n";
	outfile<<"uint64_t* forward(uint64_t* c)\n{\n\tuint64_t z[] = "<<format(z)<<";\n\t";
	outfile<<"uint64_t temp;\n\tuint64_t k=1;\n\tuint64_t j=0;\n\n";
	
	while(l>0)		//remove while for single loop testiing <-----
	{
        cout<<l<<" loop factor: ";
		cin>>factor;
		
		if(factor > n/2 || factor%2 != 0 && factor != 1)		//factor validation check
		{
			cout<<"invalid factor"<<endl;
			continue;
		}

		if(factor == n/2)	//fully unrolled
			for(int s=0; s<n; s=j+l)
			{
				for(j=s; j<s+l; j++)
				{
					outfile<<"\ttemp = (z[k] * c["<<j+l<<"])%"<<q<<";\n";
					outfile<<"\tc["<<j+l<<"] = (c["<<j<<"]+"<<q<<" - temp)%"<<q<<";\n";
					outfile<<"\tc["<<j<<"] = (c["<<j<<"] + temp)%"<<q<<";\n\n";
				}
				outfile<<"\tk++;\n\n";
			}

        else if(factor < l)	//unroll j
        {
            outfile<<"\tfor(uint64_t s=0; s<"<<n<<"; s=j+"<<l<<")\n\t{\n";
            outfile<<"\t\tfor(j=s; j<s+"<<l<<"; j+="<<factor<<")\n\t\t{\n";

            for(int i=0; i<factor; i++)
            {
                outfile<<"\t\t\ttemp = (z[k] * c[j+"<<i+l<<"])%"<<q<<";\n";
                outfile<<"\t\t\tc[j+"<<i+l<<"] = (c[j+"<<i<<"] + "<<q<<" - temp)%"<<q<<";\n";
                outfile<<"\t\t\tc[j+"<<i<<"] = (c[j+"<<i<<"] + temp)%"<<q<<";\n";         
                i == (factor-1) ? outfile<<"\t\t}\n\n\t\tk++;\n\t}\n\n":outfile<<"\n";
            }
        }   

        else		//unroll j and s
        {
            outfile<<"\tfor(uint64_t s=0; s<"<<n<<"; s+="<<2*l<<")\n\t{\n";
			for(int i=0; i<l; i++)
			{
				outfile<<"\t\ttemp = (z[k] * c[s+"<<l+i<<"])%"<<q<<";\n";
				outfile<<"\t\tc[s+"<<l+i<<"] = (c[s+"<<i<<"] + "<<q<<" - temp)%"<<q<<";\n";
				outfile<<"\t\tc[s+"<<i<<"] = (c[s+"<<i<<"] + temp)%"<<q<<";\n";
			}
			outfile<<"\t\tk++;\n\n";
			outfile<<"\t}\n\n";
        }
		
        l/=2;
	}

	outfile<<"\treturn c;\n}\n\nuint64_t c[] = "<<format(c)<<";\n\n";
	outfile<<"int main()\n{\n\tuint64_t total = 0;\n\tfor(int i=0; i<1000; i++)\n\t{\n\t\tuint64_t time = rdtscp64();\n\t\tforward(c);\n\t\ttotal += rdtscp64()-time;\n\t}\n\n\tcout<<total/1000<<endl;\n\treturn 0;\n}";
	outfile.close();
	
	return 0;
}		

//g++ *n.cpp -O3 -fno-unroll-loops -o nO.o &&  g++ *n.cpp -fno-unroll-loops -o n.o
// ./nO.o && ./n.o


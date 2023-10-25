/*
Levi Pelham - 1695061
Post-Quatum Cryptography - Optimising the Number Theoretic Transform 
Master of Cyber Security Thesis
The University of Adelaide
compiliation -> g++ *n.cpp -O3 -fno-unroll-loops -o nO.o &&  g++ *n.cpp -fno-unroll-loops -o n.o
*/

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace std;

#define maxBits 24 //maximum lazy reduction bits
#define seed 0 //set random seed

//used to calculate the execution cycles
static inline uint64_t rdtscp64() 
{
	uint32_t low, high;
	asm volatile ("rdtscp": "=a" (low), "=d" (high) :: "ecx");
	return (((uint64_t)high) << 32) | low;
}

//returns the inverse modulus
uint64_t inverse(uint64_t a, uint64_t mod)	
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

//used to precompute the roots of unity for zeta
vector<uint64_t> precompute(uint64_t n, uint64_t zeta, uint64_t q, int flag)		
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

//generate n random values within the modulus field q
vector<uint64_t> generate(uint64_t n, uint64_t q)
{
	vector<uint64_t> t;
	for(int i=0; i<n; i++)
		t.push_back(rand()%q);
	return t;
}

//testing purposes, print the vector
void print(vector<uint64_t> a)
{
	for(int i=0; i<a.size(); i++)
		cout<<a.at(i)<<" ";
	cout<<endl;
}

//formats a vector to write to file
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

//used to format the file name to be generated
string filename(vector<uint64_t> a, int n, int l)
{
	stringstream ss;
	if(a.size() == 1)
		ss<<"data/single/("<<l<<")";
	else
		ss<<"data/"<<n<<"/";

	for(int i=0; i<a.size(); i++)
	{
		ss<<a.at(i);
		i != a.size()-1 ? ss<<"," : ss<<".cpp";
	}

	return ss.str();
}

void fullUnroll(vector<uint64_t> c, uint64_t q, vector<uint64_t> z, ofstream &outfile)
{
	uint64_t k = 1;
	uint64_t j = 0;
	uint64_t n = c.size();
	uint64_t tbits;
	uint64_t initBits = ceil(log2(q));
	vector<uint64_t> bits(n,initBits);
	uint64_t bk = (uint64_t)1<<32;
	float rik = floor(bk/q)*pow(bk,-1);

	for(int l=n/2; l>0; l/=2)			//forward NTT			
		for(int s=0; s<n; s=j+l)
		{
			for(j=s; j<s+l; j++)
			{
				outfile<<"\t\ttemp = ("<<z[k]<<"* c["<<j+l<<"])%"<<q<<";\n";
				tbits = initBits;


				/*Lazy redcution default modulus*/
				// if(bits[j]+1 > maxBits)
				// {
				// 	outfile<<"\t\tc["<<j+l<<"] = (c["<<j<<"] + "<<q<<" - temp)%"<<q<<";\n";
				// 	bits[j+l] = initBits; 
				// 	outfile<<"\t\tc["<<j<<"] = (c["<<j<<"] + temp)%"<<q<<";\n\n";
				// 	bits[j] = initBits;
				// }
				// else
				// {
				// 	outfile<<"\t\tc["<<j+l<<"] = (c["<<j<<"] + "<<q<<" - temp);\n";
				// 	bits[j+l]++; 
				// 	outfile<<"\t\tc["<<j<<"] += temp;\n\n";
				// 	bits[j]++;
				// }

				/*Lazy Reduction with barrett reduction*/
				outfile<<"\ttemp = ("<<z[k]<<" * c["<<j+l<<"]) - "<<q<<" * floor(("<<z[k]<<" * c["<<j+l<<"]) * "<<rik<<");\n";
				tbits = initBits;
				
				if(bits[j]+1 > maxBits)									
				{
					outfile<<"\tc["<<j+l<<"] = (c["<<j<<"] + "<<q<<" - temp) - "<<q<<" * floor((c["<<j<<"] + "<<q<<" - temp) * "<<rik<<");\n";
					bits[j+l] = initBits; 
					outfile<<"\tc["<<j<<"] = (c["<<j<<"] + temp) - "<<q<<" * floor((c["<<j<<"] + temp) * "<<rik<<");\n\n";
					bits[j] = initBits;
				}

				else
				{
					outfile<<"\t"<<"c["<<j+l<<"] = (c["<<j<<"] + "<<q<<" - temp);\n";
					bits[j+l]++;
					outfile<<"\t"<<"c["<<j<<"] = (c["<<j<<"] + temp);\n\n";
					bits[j]++;
				}
			}

			k++;
		}	

		outfile<<"\tfor(int i=0; i<"<<n<<"; i++)\n\t\tc[i] %= "<<q<<";\n\n";
}

int main()	
{
	uint64_t q,zeta,n;		//enter parameters
	srand(0);

	n = 1024;
	q = 11777;
	zeta = 5558;

	vector<uint64_t> z = precompute(n,zeta,q,0);			//precompute all values of zeta and inverse zeta
	vector<uint64_t> invz = precompute(n,zeta,q,1);
	vector<uint64_t> c = generate(n,q);
	uint64_t inv = inverse(n,q);

	vector<uint64_t> factor = {512,512,512,512,512,512,512,512,512,512};	//unrolling factor
	uint64_t l = n/2;	//l loop
	ofstream outfile;			//generate and populate file
	outfile.open(filename(factor,n,l));
	cout<<filename(factor,n,l)<<endl;
	uint64_t k = 1;
	uint64_t j=0;
	uint64_t index = 0;

	//file structure
	outfile<<"#include <iostream>\n#include <fstream>\n#include <sstream>\n#include <cmath>\nusing namespace std;\n\n";
	outfile<<"static inline uint64_t rdtscp64() {\n\tuint32_t low, high;\n\tasm volatile (\"rdtscp\": \"=a\" (low), \"=d\" (high) :: \"ecx\");\n\treturn (((uint64_t)high) << 32) | low;\n}\n\n";
	outfile<<"uint64_t* forward(uint64_t* c)\n{\n\tuint64_t z[] = "<<format(z)<<";\n\t";
	outfile<<"uint64_t temp;\n\tuint64_t k=1;\n\tuint64_t j=0;\n\n";

	bool flag = true;

	//check if fully unrolled to use lazy reduction
	for(auto x:factor)
		if(x != n/2)
			flag = false;
		
	if(flag == true)
		fullUnroll(c,q,z,outfile);


	//fully rolled or partially unrolled
	else
		while(l>0)
		{
			if(factor[index] == n/2)	//fully unrolled
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

			else if(factor[index] < l)	//unroll j
			{
				outfile<<"\tfor(uint64_t s=0; s<"<<n<<"; s=j+"<<l<<")\n\t{\n";
				outfile<<"\t\tfor(j=s; j<s+"<<l<<"; j+="<<factor[index]<<")\n\t\t{\n";

				for(int i=0; i<factor[index]; i++)
				{
					outfile<<"\t\t\ttemp = (z[k] * c[j+"<<i+l<<"])%"<<q<<";\n";
					outfile<<"\t\t\tc[j+"<<i+l<<"] = (c[j+"<<i<<"] + "<<q<<" - temp)%"<<q<<";\n";
					outfile<<"\t\t\tc[j+"<<i<<"] = (c[j+"<<i<<"] + temp)%"<<q<<";\n";         
					i == (factor[index]-1) ? outfile<<"\t\t}\n\n\t\tk++;\n\t}\n\n":outfile<<"\n";
				}
			}   

			else		//unroll j and s
			{
				outfile<<"\tfor(uint64_t s=0; s<"<<n<<"; s+="<<2*l<<")\n\t{\n";
				
				for(int i=0; i<l; i++)
				{
					outfile<<"\t\ttemp = (z[k] * c[s+"<<l+i<<"])%"<<q<<";\n";
					outfile<<"\t\tc[s+"<<l+i<<"] = (c[s+"<<i<<"] + "<<q<<" - temp)%"<<q<<";\n";
					outfile<<"\t\tc[s+"<<i<<"] = (c[s+"<<i<<"] + temp)%"<<q<<";\n\n";
				}

				outfile<<"\t\tk++;\n";
				outfile<<"\t}\n\n";
				
			}
			
			l/=2;
			index++;
		}


	//main function and results
	outfile<<"\treturn c;\n}\n\nuint64_t c[] = "<<format(c)<<";\n\n";
	outfile<<"int main()\n{\n\tuint64_t total = 0;\n\tofstream results;\n\tresults.open(\"results/"<<format(factor)<<".txt\");\n\n\t";
	outfile<<"for(int i=0; i<200000; i++)\n\t{\n\t\tuint64_t time = rdtscp64();\n\t\tforward(c);\n\t\tresults<<rdtscp64()-time<<endl;\n\t}\n\n\tresults.close();\n\treturn 0;\n}";
	outfile.close();

	return 0;
}		



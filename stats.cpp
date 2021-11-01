#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <dirent.h>
using namespace std;

int main()
{
    DIR *dir;
    struct dirent *ent;     //open current directory and skip current and parent nodes
    dir = opendir(".");
    ent = readdir(dir);     
    
    while(ent->d_name[0] != '{')
        ent = readdir(dir);
    
    ofstream outfile;
    outfile.open("statistics.txt");
    outfile<<"file_name:\t\t\t\t\t95% CI\t Variance\t\tQ1 : Median : Q3"<<endl;
    string buffer = "";

    while(ent != NULL) //iterate over directory
    {
        ifstream infile;
        infile.open(ent->d_name);
        vector<int> times;
        long long total = 0;    //total value of data points
        int dataPoints = 0;    //number of datd points
        long long sdTemp = 0;   //standard deviation before sqrt

        while(getline(infile, buffer))
        {
            int temp = stoi(buffer);
            total+=temp;                
            dataPoints++;
            times.push_back(temp);
        }

        int mean = total/dataPoints;    //calculate mean of data

        for(int i=0; i<times.size(); i++)
            sdTemp += pow(times.at(i)-mean,2);
        
        sdTemp/=dataPoints-1;
        int sd = sqrt(sdTemp);

        int lci = mean-(1.96*sd)/sqrt(dataPoints); 
        int hci = mean+(1.96*sd)/sqrt(dataPoints); 
        outfile<<ent->d_name<<"\t"<<lci<<"\t"<<hci<<"\t("<<hci-lci<<")\t"<<times[50000]<<"\t"<<times[100000]<<"\t"<<times[150000]<<endl;
        
        infile.close();
        ent = readdir(dir);
    }
    
    outfile.close();
    closedir(dir);
 

    return 0;
}
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <dirent.h>
#include <sstream>
using namespace std;

struct node
{
    int low;
    int high;
    int Q1;
    int Q3;
    string fileName;
};

void printNode(vector<node> a)
{
	for(int i=0; i<a.size(); i++)
		cout<<a.at(i).fileName<<" : "<<a.at(i).low<<" : "<<a.at(i).high<<" : "<<a.at(i).Q1<<" : "<<a.at(i).Q3<<endl;
	cout<<endl;
}

void fillVector(vector<node>* vec)
{
    ifstream stats;
    stats.open("statistics.txt");
    string buffer = "";
    getline(stats, buffer);
    node a;

    while(getline(stats, buffer, '\t'))
    {
        a.fileName = buffer;
        getline(stats, buffer, '\t');
        a.low = stoi(buffer);
        getline(stats, buffer, '\t');
        a.high = stoi(buffer);
        getline(stats, buffer, '\t');
        getline(stats, buffer, '\t');
        a.Q1 = stoi(buffer);
        getline(stats, buffer, '\t');
        getline(stats, buffer, '\n');
        a.Q3 = stoi(buffer);        
        vec->push_back(a);
    }
    stats.close();
    return;
}

int main()
{
    vector<node> a;
    fillVector(&a);

    DIR *dir;
    struct dirent *ent;     //open current directory and skip current and parent nodes
    dir = opendir(".");
    ent = readdir(dir);     
    
    while(ent->d_name[0] != '{')
        ent = readdir(dir);

    stringstream t;
    string buffer = "";
    int index = 0;
    while(ent != NULL) //iterate over directory
    {
        ifstream infile;
        infile.open(ent->d_name);
        int count = 0;
        t <<"IQR/IQR_"<<ent->d_name;
        ofstream outfile;
        outfile.open(t.str());

        while(getline(infile, buffer))
        {
            int value = stoi(buffer);
            // if(value >= a.at(index).low && value <= a.at(index).high)
            //     outfile<<value<<endl;

            if(count >= 50000 && count <= 150000)
                outfile<<value<<endl;

            count++;
        }

        index++;
        infile.close();
        outfile.close();
        t.str("");
        ent = readdir(dir);
    }

    closedir(dir);
 
    return 0;
}
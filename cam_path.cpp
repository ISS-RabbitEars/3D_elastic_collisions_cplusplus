#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

#define PI 3.141592653589793238462643383

int main()
{
	int N;
	double r,z,theta,delta;
    ifstream in;
    ofstream file;
	
    in.open("pos/nframes.dat");
    in>>N;
    in.close();
    
    r=18;
    z=0;
    delta=2*PI/(N-1);
	theta=0;
    
    string pts="pos/cam/cam_";
	for(int i=0;i<N;i++)
	{
		theta=i*delta;
        string pu,snum;
        stringstream ssnum;
        ssnum<<i+1;
        snum=ssnum.str();
        ssnum.str("");
        pu=pts+snum;
        file.open(pu.c_str());
		file<<"<"<<r*cos(theta)<<", "<<r*sin(theta)<<", "<<z<<">"<<endl;
        file.close();
	}
    
	return 0;
}

#undef PI


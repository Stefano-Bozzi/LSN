#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;


double Error(double sum, double sum2, int iblk){
    if (iblk!=0) return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
    if (iblk==0) return 0;
}

int main(int argc, char* argv[]){

	string name=argv[1];
//lettura da file
	ifstream ifile;
	ofstream err;
	ifile.open(name);
	cout<<"aperto file: "<<name; 
	//check to see that the file was opened correctly:
	if (!ifile.is_open()) {
		cerr << "There was a problem opening the input file: "+name+"\n";
		exit(1);
	}
	double av=0.0;
	double av2=0.0;
	double sum=0.0;
	double sum2=0.0;
	int iblk=0;
	double Err=0.0;
	while (ifile >>av>>av) {
		//cout<<"lettura da file: "<<name<<endl; 
		av2=av*av;
		sum=sum+av;
		sum2=sum2+av2;
		Err=Error(sum/iblk, sum2/iblk, iblk);
		err.open("error.dat",ios::app);
		//cout<<"opened"<<"error"<<name<<endl;
		err<<Err<<"\n";
		cout<<Err<<"\n";
		err.close();
		iblk++;
	}

	
	return 0;
}


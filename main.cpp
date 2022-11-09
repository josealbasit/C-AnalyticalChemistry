#define _USE_MATH_DEFINES
#include<cmath>
#include <stdlib.h>
#include <stdio.h>
#include<iostream>
#include<Eigen/Dense>
#include<fstream>
#include<iostream>
#include <stdlib.h>     //for using the function sleep
#include<string>
#include <vector>


using namespace Eigen;
using namespace std;

void saveData(string fileName, MatrixXd matrix)
{
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ",", "\n");
	ofstream file(fileName);
	if (file.is_open())
	{
		file << matrix.format(CSVFormat);
		file.close();
	}
}

MatrixXd openData(string fileToOpen)
{
	vector<double> matrixEntries;
	// in this object we store the data from the matrix
	ifstream matrixDataFile(fileToOpen);

	// this variable is used to store the row of the matrix that contains commas 
	string matrixRowString;

	// this variable is used to store the matrix entry;
	string matrixEntry;

	// this variable is used to track the number of rows
	int matrixRowNumber = 0;


	while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
	{
		stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

		while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
		{
			matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
		}
		matrixRowNumber++; //update the column numbers
	}

	// here we convet the vector variable into the matrix and return the resulting object, 
	// note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
	if (matrixRowNumber != 0)
	{
		return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
	}
	else {
		MatrixXd a;

		return a;
	}
}

class Chromatogram {
public:
	//Attributes
	float deadTime;
	float resolution;
	MatrixXd signal;
	double phi; // Organic solvent
	double additive; // Additive concentration 
	double chrom_resolution;
    //Methods 
	double calculate_resolution(MatrixXd peak_areas) {}; //This function calculates resolution of the chromatogram based on peak-resolution 
};

class Peak {
public:
	//Attributes
	double retention; //Retention factor
	double tR; //Retention time
	double N;  //Theoretical plates
	MatrixXd signal; //Numerical signal
	double sigma; //Standard deviation of the peak
	double A; //Left half-width 
	double B; //Right half-width
	float Area;



	//Methods
	VectorXd calculate_signalSIGMA(double delta_t)
	{
		double tmax = 300;
		VectorXd t = VectorXd::LinSpaced(tmax/delta_t, 0, tmax);
		RowVectorXd peak_signal(size(t));
		for (int i = 0; i <= size(t); i++)
		{
			double calc = (1 / sigma)*(2)* exp(-0.5*pow((t(i) - tR / sigma), 2)); //Signal equation
			peak_signal(i) = calc;
		}
		return peak_signal;
	};
	VectorXd calculate_signalPMG1(double delta_t) //Model PMG1
	{
		//Parameters needed for signal declaration
		double tmax = 300; 
		const double F=sqrt(-2*log(0.1));
		double s0 = B * (1 - (B - A) / (A + B)) / F; 
		double s1 = (B - A) / ((A + B)*F);  
		double hmax = 1 / (sqrt(2 * M_PI)*(A + B)); 
		double b1 = (s0*A) / pow((s0 + s1 * (-A)), 3);
		double a1 = (hmax / 10)*exp(b1*A);
		double h = tmax / delta_t; //Signal pace
		VectorXd t1 = VectorXd::LinSpaced(h, 0, tR-A);
		VectorXd t = VectorXd::LinSpaced(h, tR - A, tR + B); 
		VectorXd t2 = VectorXd::LinSpaced(h, tR +B, tmax);
		RowVectorXd peak_signal1(size(t1)), peak_signal(size(t)), peak_signal2(size(t2));
		for (int i = 0; i <= size(t1); i++) {
			double calc= a1*exp(b1*(t1(i) - tR));
			peak_signal1(i) = calc; 
		};
		for (int i = 0; i <= size(t); i++) {
			double calc = a1 * exp(b1*(t(i) - tR));
			peak_signal(i) = calc;
		};
		for (int i = 0; i <= size(t); i++) {
			double calc = a1 * exp(b1*(t(i) - tR));
			peak_signal(i) = calc;
		};

	    
	};

public:
	float resolutionCalculation() {

	};
};






int main() {
	string fileNameRetention;
	string fileNameA;
	string fileNameB;
	string fileNameParameters;
	Chromatogram A;
	cout << "Introduce dead time associated with this chromatogram"<< endl;
	cin >> A.deadTime;
	cout << A.deadTime << endl;
	cout << "Intro retention data file name" << endl;
	cin >> fileNameRetention;
	cout << "Intro A10 data file name" << endl;
	cin >> fileNameA;
	cout << "Intro B10 data file name" << endl;
	cin >> fileNameB;
	cout << "Intro Parameters data file name" << endl;
	cin >> fileNameParameters; 
	cout << "/n Intro coefficients of your equation /n " << endl; // At this moment we only support one type of equation
	MatrixXd RetentionData = openData(fileNameRetention);
	MatrixXd AData = openData(fileNameA);
	MatrixXd BData = openData(fileNameB);
	MatrixXd ParameterData = openData(fileNameParameters);
	cout << RetentionData << endl;
	cout << AData << endl; 
	cout << BData << endl;
	cout << ParameterData << endl; 
	int m = RetentionData.rows();
	int n = int(RetentionData.cols() - 1);
	cout << RetentionData.size() << n << endl;
	system("PAUSE");
	//MatrixXd RetentionTime = (DeadTime.replicate(1,n)).cwiseProduct(Retention) + DeadTime.replicate(1,n);
	//cout << RetentionTime << endl; 

	return 0;

}
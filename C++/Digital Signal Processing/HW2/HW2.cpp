////////////////////////////////////
// ECE529 HW2 C++ Problem         //
//                                //
// John Claus 02/14/2017          //
////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>

using namespace std;

//Global Variables - Global Arrays used in this program because pointers are required in C++ for returning arrays from functions
float x[5000], y[5000], a[5000], b[5000], h[5000];
int stag = 0;
ofstream myfile;



//Functions
void filter();
void conv();
void InitY();
void YtoH();
void DisplayArray(int L, char T);
void BuildArray(int L, char T);


int main()
{
	char f = 'n', s = 'n';
	string filename;
	while (f != 'Q')
	{
		cout << "(F)ilter, (C)onvolution, or (Q)uit: ";
		cin >> f;
		cout << "\n";
		if (f != 'Q')
		{
			cout << "Save to  File?(Y/N): ";
			cin >> s;
			cout << "\n";
		}
		if (s == 'Y')
		{
			stag = 1;
			cout << "Enter Filename: ";
			cin >> filename;
			cout << "\n";
			myfile.open(filename);
			if (f == 'F') filter();
			if (f == 'C') conv();
			myfile.close();
		}
		if (s != 'Y')
		{
			stag = 0;
			if (f == 'F') filter();
			if (f == 'C') conv();
		}
	}

}


//Function is used to build arrays to be used in other functions
void BuildArray(int L, char T)
{
	float val;
	for (int i = 0; i <= L; i++)
	{
		if ((T != 'a') || ((T == 'a') && (i != 0)))
		{
			cout << "Enter " << T << "[" << i << "]: ";
			cin >> val;
			if (T == 'x') x[i] = val;
			if (T == 'y') y[i] = val;
			if (T == 'a') a[i] = val;
			if (T == 'b') b[i] = val;
			if (T == 'h') h[i] = val;
		}
	}
	cout << "\n";
	return;
}

//Function is used to display array values 
void DisplayArray(int L, char T)
{
	if (stag == 1) myfile << T << "[n] = ";
	cout << T << "[n] = ";
	for (int i = 0; i <= L; i++)
	{
		if (T == 'x')
		{
			cout << x[i] << " ";
			if (stag == 1) myfile << x[i] << " ";
		}
		if (T == 'y')
		{
			cout << y[i] << " ";
			if (stag == 1) myfile << y[i] << " ";
		}
		if (T == 'a')
		{
			cout << a[i] << " ";
			if (stag == 1) myfile << a[i] << " ";
		}
		if (T == 'b')
		{
			cout << b[i] << " ";
			if (stag == 1) myfile << b[i] << " ";
		}
		if (T == 'h')
		{
			cout << h[i] << " ";
			if (stag == 1) myfile << h[i] << " ";
		}
	}
	cout << "\n \n";
	if (stag == 1) myfile << "\n \n";
	return;
}

//Function is used to initialize the output array to zero
void InitY()
{
	for (int n = 0; n < 5000; n++)
		y[n] = 0;
}

//Function is used to make the output array a response
void YtoH()
{
	for (int n = 0; n < 5000; n++)
		h[n] = y[n];
}

/*Function is used to perform filtering as described in HW#2 Part A:

This function implements an arbitrary LCCDE according to y[n] = -SUM_{ k = 1 }^{N} a_k y[n - k] + SUM_{ k = 0 }^{M} b_k x[n - k]
The input signal x[n] is given for 0 <= n < L.Assume zero initial conditions.The main program allocates the four arrays and fills the
b, a, and x arrays.This function then computes y[n] for 0 <= n <L.
Inputs: b[], M, a[], N, x[], L.
Output : y[]

*/


void filter()
{
	int n, k, L, N, M;
	float xtemp, ytemp;
	char c;

	//The max index is requested from the user
	cout << "Enter a desired L where (0 <= n < L): ";
	cin >> L;
	cout << "\n";

	//This section is where the output array y[n] is initialized and x[n] is created for the index size L
	InitY();
	DisplayArray((L-1), 'y');
	BuildArray((L-1), 'x');
	DisplayArray((L-1), 'x');
	cout << "\n";

	//This is where difference equation setting are requested from the user
	cout << "Enter largest X offset: ";
	cin >> M;
	cout << "\n";
	cout << "Enter largest Y offset: ";
	cin >> N;
	cout << "\n";


	//This section is where the constant arrays a and b are created for the index size L
	BuildArray(N, 'a');
	//DisplayArray(N, 'a');
	BuildArray(M, 'b');
	//DisplayArray(L, 'b');


	//This section returns the difference equations as they would appear for values M and N
	cout << "\n\n" << "y[n] = ";
	if (stag == 1) myfile << "\n" << "y[n] = ";
	if (M >= 0)
	{
		cout << b[0] <<"*x[n]";
		if (stag == 1) myfile << b[0] << "*x[n]";
		if (M > 0)
			for (int i = 1; i <= M; i++)
			{
				cout << " + " << b[i] << "*x[n-" << i << "]";
				if (stag == 1) myfile << " + " << b[i] << "*x[n-" << i << "]";
			}
	}
	if (N >= 0)
	{
		if (N > 0)
			for (int i = 1; i <= N; i++)
			{
				cout << " - " << a[i] << "*y[n-" << i << "]";
				if (stag == 1) myfile << " - " << a[i] << "*y[n-" << i << "]";
			}
	}
	cout << "\n";
	if (stag == 1) myfile << "\n";


	//This section is where the filtering is performed
	for (n = 0; n < L; n++)
	{
		xtemp = 0;
		//Here the x components from k=0 to M are summed for each n
		for (k = 0; k <= M; k++)
		{
			if ((n - k) >= 0) xtemp = xtemp + (b[k]*x[(n - k)]);
			//Below was used for debugging
			/*cout << "xtemp = " << xtemp;
			cout << ", x[n-k] = " << x[(n - k)];
			cout << ", n = " << n;
			cout << ", k = " << k << "\n";*/
		}
		//cout << "\n";
		ytemp = 0;
		//Here the y components from k=0 to N are summed for each n
		for (k = 0; k <= N; k++)
		{
			if ((n - k) >= 0) ytemp = ytemp + (a[k]*y[(n - k)]);
			//Below s used for debugging
			/*cout << "ytemp = " << ytemp;
			cout << ", y[n-k] = " << x[(n - k)];
			cout << ", n = " << n;
			cout << ", k = " << k << "\n";*/
		}
		//Here those summed values are combined for each y[n]
		y[n] = xtemp - ytemp;
		//cout << "xtemp = " << xtemp << ", ytemp =" << ytemp << ", y[n] = " << y[n] << "\n";  //more dubugging
	}

	//This section renames y[n] as h[n] if desired
	cout << "Make y[n] be h[n]? (Y/N): ";
	cin >> c;
	if (c == 'Y')
	{
		YtoH();
		cout << "\n" << "Solution:\n";
		if (stag == 1) myfile << "\n" << "Solution:\n";
		DisplayArray(L, 'h');
	}
	if (c != 'Y')
	{
		cout << "\n" << "Solution:\n";
		if (stag == 1) myfile << "\n" << "Solution:\n";
		DisplayArray(L, 'y');
	}
	return;
}


void conv()
{
	char r;
	int L;
	int M;
	cout << "Enter the length of the input signal L: ";
	cin >> L;
	cout << "\n";
	InitY();
	DisplayArray((L - 1), 'y');
	BuildArray((L - 1), 'x');
	DisplayArray((L - 1), 'x');
	cout << "\n";
	cout << "Have you created a response h(x)? (Y/N) ";
	cin >> r;
	cout << "\n";
	cout << "Enter the length of h[x], M: ";
	cin >> M;
	cout << "\n";
	if (r != 'Y')
	{
		YtoH();
		BuildArray(M, 'h');
	}
	DisplayArray(M, 'h');
	cout << "\n\n";
	if (stag == 1) myfile << "\n\n";
	for (int n = 0; n < L; n++)
	{
		for (int m = 0; m <= M; m++)
		{
			y[n] = y[n] + (h[n] * x[(n - m)]);
		}
	}
	cout << "Solution: \n\n";
	if (stag == 1) myfile << "Solution: \n\n";
	DisplayArray((L + M - 1), 'y');
}
#include "stdafx.h"
#include "Spline.h"
#include "utils.h"
#include <iostream>

void splineknots::Spline::Print()
{
	using namespace std;
	cout << "---------- Knot matrix ----------" << endl;
	for (size_t i = 0; i < RowsCount(); i++)
	{
		cout << "Row " << i << " :\n";
		for (size_t j = 0; j < ColumnsCount(); j++)
		{
			cout << j << ":\n"
				<< "z: " << z[i][j] << '\n'
				<< "dx: " << dx[j][i] << '\n'
				<< "dy: " << dy[i][j] << '\n'
				<< "dxy: " << dxy[i][j] << '\n';
		}
		cout << endl;
	}
	cout << "-------------------------------" << endl;
}


splineknots::Spline splineknots::Spline::EmptySpline()
{
	Spline nullval;
	return nullval;
}

bool splineknots::Spline::IsEmpty()
{
	if (RowsCount() < 1 || ColumnsCount() < 1)
		return true;
	return false;
}

splineknots::Spline::Spline()
	: x(), y(), z(), dx(), dy(), dxy()
{
}

splineknots::Spline::Spline(KnotVector rowVector, KnotVector columnVector)
	: x(std::move(rowVector)), y(std::move(columnVector)), z(), dx(), dy(), dxy()
{
	z.resize(x.size());
	dx.resize(x.size());
	dy.resize(x.size());
	dxy.resize(x.size());

	for (int i = 0; i < x.size(); ++i) {
		z[i].resize(y.size());
		dx[i].resize(y.size());
		dy[i].resize(y.size());
		dxy[i].resize(y.size());
	}
}


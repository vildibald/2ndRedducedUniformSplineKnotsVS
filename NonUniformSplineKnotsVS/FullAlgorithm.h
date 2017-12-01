//
// Created by Viliam on 15.2.2017.
//
#include "MathFunction.h"
#include "Spline.h"
#include <functional>
#include <vector>
#include "Tridiagonal.h"
#include <omp.h>
#include "Timer.h"


class FullAlgorithm
{
	InterpolativeMathFunction function_;
	Tridiagonals xTridiagonals_;
	Tridiagonals yTridiagonals_;
	Timer timer_;
	bool isParallel_;

public:
	FullAlgorithm(const InterpolativeMathFunction f);

	Spline Calculate(const KnotVector xVector,
	                 const KnotVector yVector);

	double ExecutionTime();

	double AllTime();

	void InParallel(bool value);

	bool IsParallel() const;

private:
	void Initialize(Spline& values);

	void FillDx(Spline& values);

	void FillDy(Spline& values);

	void FillDxy(Spline& values);

	void FillDyx(Spline& values);

	void Parallelize(bool inParallel);

	void InitializeTridiagonals(Spline& spline);
};

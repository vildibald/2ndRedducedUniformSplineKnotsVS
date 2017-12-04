#pragma once
#include "MathFunction.h"
#include "Tridiagonal.h"

class SplineAlgorithm
{
protected:
	InterpolativeMathFunction function_;
	Tridiagonals xTridiagonals_;
	Tridiagonals yTridiagonals_;
	Timer timer_;
	bool isParallel_;

public:

	virtual ~SplineAlgorithm()
	{
	}

	SplineAlgorithm(const SplineAlgorithm& other) = delete;
	SplineAlgorithm(SplineAlgorithm&& other) noexcept = delete;
	SplineAlgorithm& operator=(const SplineAlgorithm& other) = delete;
	SplineAlgorithm& operator=(SplineAlgorithm&& other) noexcept = delete;


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

	virtual void InitializeTridiagonals(Spline& spline) = 0;

	template <typename ParameterGetter, typename DerivationGetter, typename DerivationSetter>
	void FillD(const size_t systemCount, const size_t derivationCount, Tridiagonals& tridiagonals,
		const ParameterGetter& p, const DerivationGetter& dget, DerivationSetter& dset)
	{
		
	}
};


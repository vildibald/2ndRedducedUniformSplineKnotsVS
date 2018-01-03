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
#include "utils.h"

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

	template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
	          DerivationSetter>
	void FillD(const size_t systemCount, const size_t derivationCount, Tridiagonals& tridiagonals,
	           const DifferenceGetter h, const ParameterGetter& p, const DerivationGetter& dget,
	           DerivationSetter& dset)
	{
		utils::For(0, static_cast<int>(systemCount), 1, isParallel_, [&](int j)
		{
			Solve(systemCount, derivationCount, tridiagonals, h, p, dget, dset, j);
		});
	}

	template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
	          DerivationSetter>
	void Solve(const size_t systemCount, const size_t derivationCount, Tridiagonals& tridiagonals,
	           const DifferenceGetter h, const ParameterGetter& p, const DerivationGetter& dget,
	           DerivationSetter& dset, size_t systemIdx)
	{
		auto& tridiagonal = tridiagonals.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto threeDivH = 3.0 / tridiagonal.H();
		for (size_t i = 0; i < rightSide.size(); ++i)
		{
			rightSide[i] = threeDivH * (p(i + 2, systemIdx) - p(i, systemIdx));
		}

		rightSide.front() -= dget(0, systemIdx);;
		rightSide.back() -= dget(derivationCount - 1, systemIdx);
		tridiagonal.Solve();

		for (auto i = 0; i < rightSide.size(); ++i)
		{
			dset(i + 1, systemIdx, rightSide[i]);
		}
	}
};

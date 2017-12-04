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
		const auto& r2 = tridiagonal.RhsCoeficients()[2];
		const auto& r1 = tridiagonal.RhsCoeficients()[1];
		const auto& r0 = tridiagonal.RhsCoeficients()[0];
		const auto l2 = tridiagonal.Lhs2Coeficients().front();
		const auto l0 = tridiagonal.Lhs0Coeficients().back();

		for (size_t i = 0; i < rightSide.size(); ++i)
		{
			rightSide[i] = r2[i] * p(i + 2, systemIdx) + r1[i] * p(i + 1, systemIdx) + r0[i] * p(i, systemIdx);
		}

		rightSide.front() -= l0 * dget(0, systemIdx);;
		rightSide.back() -= l2 * dget(derivationCount - 1, systemIdx);
		tridiagonal.Solve();

		for (auto i = 0; i < rightSide.size(); ++i)
		{
			dset(i + 1, systemIdx, rightSide[i]);
		}
	}
};

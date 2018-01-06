#pragma once
#include "Tridiagonal.h"

class FullSolver
{
public:
	FullSolver();

	template <typename ParameterGetter, typename DerivationGetter, typename
	          DerivationSetter>
	void Solve(const size_t systemCount, const size_t derivationCount, Tridiagonals& tridiagonals,
	           const ParameterGetter& p, const DerivationGetter& dget,
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

	Tridiagonal CreateTridiagonal(const KnotVector& knots);
};

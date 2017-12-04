#pragma once
#include "MathFunction.h"
#include "Spline.h"
#include <vector>
#include "Tridiagonal.h"

class ReducedAlgorithm
{
	InterpolativeMathFunction function_;
	Tridiagonals xTridiagonals_;
	Tridiagonals yTridiagonals_;
	Timer timer_;
	bool isParallel_;

public:
	ReducedAlgorithm(const InterpolativeMathFunction f);

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
		const auto unknownsCount = derivationCount;
		const auto even = unknownsCount % 2 == 0;
		const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;

		auto& tridiagonal = tridiagonals.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto& r4 = tridiagonal.RhsCoeficients()[4];
		const auto& r3 = tridiagonal.RhsCoeficients()[3];
		const auto& r2 = tridiagonal.RhsCoeficients()[2];
		const auto& r1 = tridiagonal.RhsCoeficients()[1];
		const auto& r0 = tridiagonal.RhsCoeficients()[0];
		const auto l4 = tridiagonal.Lhs2Coeficients().front();
		const auto l0 = tridiagonal.Lhs0Coeficients().back();

		for (size_t i = 0; i < equationsCount - 1; ++i)
		{
			auto i21 = 2 * (i + 1);
			rightSide[i] = r4[i] * p(i21 + 2, systemIdx) + r3[i] * p(i21 + 1, systemIdx) + r2[i] * p(i21, systemIdx) + r1[i] *
				p(i21 - 1, systemIdx) + r0[i] * p(i21 - 2, systemIdx);
		}
		rightSide.front() -= l4 * dget(0, systemIdx);
		{
			const int i = equationsCount - 1;
			const auto i21 = 2 * (i + 1);
			if (even)
			{
				// TODO: Needs support for even count.
				rightSide.back() = 0;
			}
			else
			{
				rightSide.back() = r4[i] * p(i21 + 2, systemIdx) + r3[i] * p(i21 + 1, systemIdx) + r2[i] * p(i21, systemIdx) + r1[i] *
					p(i21 - 1, systemIdx) + r0[i] * p(i21 - 2, systemIdx);
			}
			rightSide.back() -= l0 * dget(derivationCount - 1, systemIdx);
		}
		tridiagonal.Solve();

		const auto& rst2 = tridiagonal.RhsCoeficients()[7];
		const auto& rst1 = tridiagonal.RhsCoeficients()[6];
		const auto& rst0 = tridiagonal.RhsCoeficients()[5];
		const auto& rstd2 = tridiagonal.RhsCoeficients()[9];
		const auto& rstd0 = tridiagonal.RhsCoeficients()[8];
		size_t i = 0;
		for (; i < equationsCount; ++i)
		{
			const auto evenI = 2 * (i + 1);
			dset(evenI, systemIdx, rightSide[i]);

			const auto oddI = 2 * i + 1;
			dset(oddI, systemIdx,
				rst2[i] * p(oddI + 1, systemIdx) + rst1[i] * p(oddI, systemIdx) + rst0[i] * p(oddI - 1, systemIdx) +
				rstd2[i] * dget(oddI + 1, systemIdx) + rstd0[i] * dget(oddI - 1, systemIdx)
			);
		}
		const auto oddI = 2 * i + 1;
		dset(oddI, systemIdx,
			rst2[i] * p(oddI + 1, systemIdx) + rst1[i] * p(oddI, systemIdx) + rst0[i] * p(oddI - 1, systemIdx) +
			rstd2[i] * dget(oddI + 1, systemIdx) + rstd0[i] * dget(oddI - 1, systemIdx)
		);
	}
};

#pragma once
#include "Tridiagonal.h"

class SecondReducedSolver
{
public:
	SecondReducedSolver();

	template <typename ParameterGetter, typename DerivationGetter, typename
	          DerivationSetter>
	void Solve(const size_t systemCount, const size_t derivationCount, Tridiagonals& tridiagonals,
	           const ParameterGetter& p, const DerivationGetter& dget,
	           DerivationSetter& dset, size_t systemIdx)
	{
		auto& tridiagonal = tridiagonals.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto unknownsCount = derivationCount;
		const auto equationsCount = tridiagonal.EquationCount();
		const auto h = tridiagonal.H();
		const auto threeDivH = 3.0 / h;
		for (size_t i = 0; i < equationsCount - 1; ++i)
		{
			const auto i41 = 4 * (i + 1);
			rightSide[i] = threeDivH * (p(i41 + 4, systemIdx) - p(i41 - 4, systemIdx)
				- 4 * (p(i41 + 3, systemIdx) - p(i41 - 3, systemIdx))
				+ 14 * (p(i41 + 2, systemIdx) - p(i41 - 2, systemIdx))
				- 52 * (p(i41 + 1, systemIdx) - p(i41 - 1, systemIdx))
			);
		}
		const auto d0 = dget(0, systemIdx);
		rightSide.front() -= d0;
		const auto dN = dget(derivationCount - 1, systemIdx);
		{
			const size_t i = equationsCount - 1;
			const auto i41 = 4 * (i + 1);
			const auto remainder = unknownsCount % 4;
			switch (remainder)
			{
			case 3:
				break;
			case 2:
				break;
			case 1:
				rightSide.back() = threeDivH * (p(i41 + 4, systemIdx) - p(i41 - 4, systemIdx)
					- 4 * (p(i41 + 3, systemIdx) - p(i41 - 3, systemIdx))
					+ 14 * (p(i41 + 2, systemIdx) - p(i41 - 2, systemIdx))
					- 52 * (p(i41 + 1, systemIdx) - p(i41 - 1, systemIdx))
				);
				break;
			default:

				break;
			}
			rightSide.back() -= dN;
		}
		tridiagonal.Solve();

		const auto minusOneDiv14 = -1.0 / 14;
		const auto twelveDivH = 4.0 * threeDivH;
		const auto oneDiv4 = 1.0 / 4;
		//const auto threeDiv4H = oneDiv4 * threeDivH;
		for (size_t i = 0; i < equationsCount; ++i)
		{
			auto idx = 4 * (i + 1);
			dset(idx, systemIdx, rightSide[i]);

			idx -= 2;
			dset(idx, systemIdx,
			     minusOneDiv14 * (
				     threeDivH * (p(idx + 2, systemIdx) - p(idx - 2, systemIdx))
				     - twelveDivH * (p(idx + 1, systemIdx) - p(idx - 1, systemIdx))
				     - dget(idx - 2, systemIdx) - dget(idx + 2, systemIdx)
			     )
			);

			--idx;
			dset(idx, systemIdx,
				oneDiv4 * (
					threeDivH * (p(idx + 1, systemIdx) - p(idx - 1, systemIdx))
					- dget(idx - 1, systemIdx) - dget(idx + 1, systemIdx)
					)
			);

			idx += 2;
			dset(idx, systemIdx,
			     oneDiv4 * (
				     threeDivH * (p(idx + 1, systemIdx) - p(idx - 1, systemIdx))
				     - dget(idx - 1, systemIdx) - dget(idx + 1, systemIdx)
			    )
				/*threeDivH * (p(idx + 1, systemIdx) - p(idx - 1, systemIdx))
				- oneDiv4 * (dget(idx - 1, systemIdx) + dget(idx + 1, systemIdx))*/
			);
		}
		{
			const size_t i = equationsCount;
			auto idx = 4 * (i + 1);

			idx -= 2;
			dset(idx, systemIdx,
			     minusOneDiv14 * (
				     threeDivH * (p(idx + 2, systemIdx) - p(idx - 2, systemIdx))
				     - twelveDivH * (p(idx + 1, systemIdx) - p(idx - 1, systemIdx))
				     - dget(idx - 2, systemIdx) - dget(idx + 2, systemIdx)
			     )
			);

			--idx;
			dset(idx, systemIdx,
			     oneDiv4 * (
				     threeDivH * (p(idx + 1, systemIdx) - p(idx - 1, systemIdx))
				     - dget(idx - 1, systemIdx) - dget(idx + 1, systemIdx)
			     )
			);

			idx += 2;
			dset(idx, systemIdx,
			     oneDiv4 * (
				     threeDivH * (p(idx + 1, systemIdx) - p(idx - 1, systemIdx))
				     - dget(idx - 1, systemIdx) - dget(idx + 1, systemIdx)
			     )
			);
		}
	}

	Tridiagonal CreateTridiagonal(const KnotVector& knots);
};

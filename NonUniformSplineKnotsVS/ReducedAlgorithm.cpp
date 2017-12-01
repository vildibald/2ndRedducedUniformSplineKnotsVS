#include "stdafx.h"
#include "utils.h"
#include "ReducedAlgorithm.h"

ReducedAlgorithm::ReducedAlgorithm(const InterpolativeMathFunction f)
	: function_(f), xTridiagonals_(),
	  yTridiagonals_(), isParallel_(false)
{
}

Spline ReducedAlgorithm::Calculate(const KnotVector xVector,
                                   const KnotVector yVector)
{
	Spline spline{std::move(xVector), std::move(yVector)};
	Initialize(spline);
	timer_.Reset();
	timer_.Start();
	FillDx(spline);
	FillDy(spline);
	FillDxy(spline);
	FillDyx(spline);
	timer_.Stop();

	return spline;
}

void ReducedAlgorithm::FillDx(Spline& spline)
{
	const auto unknownsCount = spline.RowsCount();
	const auto even = unknownsCount % 2 == 0;
	const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;
	utils::For(0, static_cast<int>(spline.ColumnsCount()), 1, isParallel_, [&](int j)
	{
		auto& tridiagonal = xTridiagonals_.Get();
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
			rightSide[i] = r4[i] * spline.Z(i21 + 2, j) + r3[i] * spline.Z(i21 + 1, j) + r2[i] * spline.Z(i21, j) + r1[i] *
				spline.Z(i21 - 1, j) + r0[i] * spline.Z(i21 - 2, j);
		}
		rightSide.front() -= l4 * spline.Dx(0, j);
		rightSide.back() -= l0 * spline.Dx(spline.RowsCount() - 1, j);
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
				rightSide.back() = r4[i] * spline.Z(i21 + 2, j) + r3[i] * spline.Z(i21 + 1, j) + r2[i] * spline.Z(i21, j) + r1[i] *
					spline.Z(i21 - 1, j) + r0[i] * spline.Z(i21 - 2, j);
			}
		}
		tridiagonal.Solve();
		for (size_t i = 0; i < rightSide.size(); ++i)
		{
			const auto i21 = 2 * (i + 1);
			spline.SetDx(i21, j, rightSide[i]);
		}

		const auto& rst2 = tridiagonal.RhsCoeficients()[7];
		const auto& rst1 = tridiagonal.RhsCoeficients()[6];
		const auto& rst0 = tridiagonal.RhsCoeficients()[5];
		const auto& rstd2 = tridiagonal.RhsCoeficients()[9];
		const auto& rstd0 = tridiagonal.RhsCoeficients()[8];

		for (size_t i = 1, idx = 0; i < spline.RowsCount(); i += 2, ++idx)
		{
			spline.SetDx(i, j,
			             rst2[idx] * spline.Z(i + 1, j) + rst1[idx] * spline.Z(i, j) + rst0[idx] * spline.Z(i - 1, j) +
			             rstd2[idx] * spline.Dx(i + 1, j) + rstd0[idx] * spline.Dx(i - 1, j)
			);
		}
	});
}

void ReducedAlgorithm::FillDy(Spline& spline)
{
	const auto unknownsCount = spline.ColumnsCount();
	const auto even = unknownsCount % 2 == 0;
	const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;
	utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i)
	{
		auto& tridiagonal = yTridiagonals_.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto& r4 = tridiagonal.RhsCoeficients()[4];
		const auto& r3 = tridiagonal.RhsCoeficients()[3];
		const auto& r2 = tridiagonal.RhsCoeficients()[2];
		const auto& r1 = tridiagonal.RhsCoeficients()[1];
		const auto& r0 = tridiagonal.RhsCoeficients()[0];
		const auto l4 = tridiagonal.Lhs2Coeficients().front();
		const auto l0 = tridiagonal.Lhs0Coeficients().back();

		for (size_t j = 0; j < equationsCount - 1; ++j)
		{
			auto j21 = 2 * (j + 1);
			rightSide[j] = r4[j] * spline.Z(i, j21 + 2) + r3[j] * spline.Z(i, j21 + 1) + r2[j] * spline.Z(i, j21) + r1[j] *
				spline.Z(i, j21 - 1) + r0[j] * spline.Z(i, j21 - 2);
		}
		rightSide.front() -= l4 * spline.Dy(i, 0);
		rightSide.back() -= l0 * spline.Dy(i, spline.ColumnsCount() - 1);
		{
			const int j = equationsCount - 1;
			const auto j21 = 2 * (j + 1);
			if (even)
			{
				// TODO: Needs support for even count.
				rightSide.back() = 0;
			}
			else
			{
				rightSide[j] = r4[j] * spline.Z(i, j21 + 2) + r3[j] * spline.Z(i, j21 + 1) + r2[j] * spline.Z(i, j21) + r1[j] *
					spline.Z(i, j21 - 1) + r0[j] * spline.Z(i, j21 - 2);
			}
		}
		tridiagonal.Solve();
		for (size_t j = 0; j < rightSide.size(); ++j)
		{
			const auto j21 = 2 * (j + 1);
			spline.SetDy(i, j21, rightSide[j]);
		}

		const auto& rst2 = tridiagonal.RhsCoeficients()[7];
		const auto& rst1 = tridiagonal.RhsCoeficients()[6];
		const auto& rst0 = tridiagonal.RhsCoeficients()[5];
		const auto& rstd2 = tridiagonal.RhsCoeficients()[9];
		const auto& rstd0 = tridiagonal.RhsCoeficients()[8];

		for (size_t j = 1, jdx = 0; j < spline.RowsCount(); j += 2, ++jdx)
		{
			spline.SetDy(j, j,
			             rst2[jdx] * spline.Z(i, j + 1) + rst1[jdx] * spline.Z(i, j) + rst0[jdx] * spline.Z(i, j - 1) +
			             rstd2[jdx] * spline.Dy(i, j + 1) + rstd0[jdx] * spline.Dy(i, j - 1)
			);
		}
	});
}

void ReducedAlgorithm::FillDxy(Spline& spline)
{
	const auto unknownsCount = spline.RowsCount();
	const auto even = unknownsCount % 2 == 0;
	const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;
	size_t loop[] = { 0, spline.ColumnsCount() - 1 };
	for(auto j : loop)
	{
		auto& tridiagonal = xTridiagonals_.Get();
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
			const auto i21 = 2 * (i + 1);
			rightSide[i] = r4[i] * spline.Dy(i21 + 2, j) + r3[i] * spline.Dy(i21 + 1, j) + r2[i] * spline.Dy(i21, j) + r1[i] *
				spline.Dy(i21 - 1, j) + r0[i] * spline.Dy(i21 - 2, j);
		}
		rightSide.front() -= l4 * spline.Dxy(0, j);
		rightSide.back() -= l0 * spline.Dxy(spline.RowsCount() - 1, j);
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
				rightSide[i] = r4[i] * spline.Dy(i21 + 2, j) + r3[i] * spline.Dy(i21 + 1, j) + r2[i] * spline.Dy(i21, j) + r1[i] *
					spline.Dy(i21 - 1, j) + r0[i] * spline.Dy(i21 - 2, j);
			}
		}
		tridiagonal.Solve();
		for (size_t i = 0; i < rightSide.size(); ++i)
		{
			const auto i21 = 2 * (i + 1);
			spline.SetDxy(i21, j, rightSide[i]);
		}

		const auto& rst2 = tridiagonal.RhsCoeficients()[7];
		const auto& rst1 = tridiagonal.RhsCoeficients()[6];
		const auto& rst0 = tridiagonal.RhsCoeficients()[5];
		const auto& rstd2 = tridiagonal.RhsCoeficients()[9];
		const auto& rstd0 = tridiagonal.RhsCoeficients()[8];

		for (size_t i = 1, idx = 0; i < spline.RowsCount(); i += 2, ++idx)
		{
			spline.SetDxy(i, j,
				rst2[idx] * spline.Dy(i + 1, j) + rst1[idx] * spline.Dy(i, j) + rst0[idx] * spline.Dy(i - 1, j) +
				rstd2[idx] * spline.Dxy(i + 1, j) + rstd0[idx] * spline.Dxy(i - 1, j)
			);
		}
	}
}

void ReducedAlgorithm::FillDyx(Spline& spline)
{
	const auto unknownsCount = spline.ColumnsCount();
	const auto even = unknownsCount % 2 == 0;
	const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;
	utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i)
	{
		auto& tridiagonal = yTridiagonals_.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto& r4 = tridiagonal.RhsCoeficients()[4];
		const auto& r3 = tridiagonal.RhsCoeficients()[3];
		const auto& r2 = tridiagonal.RhsCoeficients()[2];
		const auto& r1 = tridiagonal.RhsCoeficients()[1];
		const auto& r0 = tridiagonal.RhsCoeficients()[0];
		const auto l4 = tridiagonal.Lhs2Coeficients().front();
		const auto l0 = tridiagonal.Lhs0Coeficients().back();

		for (size_t j = 0; j < equationsCount - 1; ++j)
		{
			const auto j21 = 2 * (j + 1);
			rightSide[j] = r4[j] * spline.Dx(i, j21 + 2) + r3[j] * spline.Dx(i, j21 + 1) + r2[j] * spline.Dx(i, j21) + r1[j] *
				spline.Dx(i, j21 - 1) + r0[j] * spline.Dx(i, j21 - 2);
		}
		rightSide.front() -= l4 * spline.Dxy(i, 0);
		rightSide.back() -= l0 * spline.Dxy(i, spline.ColumnsCount() - 1);
		{
			const int j = equationsCount - 1;
			const auto j21 = 2 * (j + 1);
			if (even)
			{
				// TODO: Needs support for even count.
				rightSide.back() = 0;
			}
			else
			{
				rightSide[j] = r4[j] * spline.Dx(i, j21 + 2) + r3[j] * spline.Dx(i, j21 + 1) + r2[j] * spline.Dx(i, j21) + r1[j] *
					spline.Dx(i, j21 - 1) + r0[j] * spline.Dx(i, j21 - 2);
			}
		}
		tridiagonal.Solve();
		for (size_t j = 0; j < rightSide.size(); ++j)
		{
			auto j21 = 2 * (j + 1);
			spline.SetDxy(i, j21, rightSide[j]);
		}

		const auto& rst2 = tridiagonal.RhsCoeficients()[7];
		const auto& rst1 = tridiagonal.RhsCoeficients()[6];
		const auto& rst0 = tridiagonal.RhsCoeficients()[5];
		const auto& rstd2 = tridiagonal.RhsCoeficients()[9];
		const auto& rstd0 = tridiagonal.RhsCoeficients()[8];

		for (size_t j = 1, jdx = 0; j < spline.RowsCount(); j += 2, ++jdx)
		{
			spline.SetDxy(j, j,
				rst2[jdx] * spline.Dx(i, j + 1) + rst1[jdx] * spline.Dx(i, j) + rst0[jdx] * spline.Dx(i, j - 1) +
				rstd2[jdx] * spline.Dxy(i, j + 1) + rstd0[jdx] * spline.Dxy(i, j - 1)
			);
		}
	});
}

void ReducedAlgorithm::Initialize(Spline& spline)
{
	InitializeTridiagonals(spline);
	spline.Initialize(function_);
}

bool ReducedAlgorithm::IsParallel() const
{
	return isParallel_;
}

void ReducedAlgorithm::InParallel(const bool value)
{
	isParallel_ = value;
}

void ReducedAlgorithm::Parallelize(const bool inParallel)
{
	xTridiagonals_.Parallelize(inParallel);
	yTridiagonals_.Parallelize(inParallel);
}

void ReducedAlgorithm::InitializeTridiagonals(Spline& spline)
{
	xTridiagonals_.GetAll().clear();
	xTridiagonals_.GetAll().clear();
	xTridiagonals_.GetAll().emplace_back(
		Tridiagonal::Factory::CreateReducedTridiagonal(spline.X, spline.RowsCount()));
	yTridiagonals_.GetAll().emplace_back(
		Tridiagonal::Factory::CreateReducedTridiagonal(spline.Y, spline.ColumnsCount()));
	Parallelize(isParallel_);
}

double ReducedAlgorithm::ExecutionTime()
{
	return timer_.ExecutionTime();
}

double ReducedAlgorithm::AllTime()
{
	return timer_.AllTime() +
		xTridiagonals_.GetAll()[0].AllTime() + yTridiagonals_.GetAll()[0].AllTime();
}

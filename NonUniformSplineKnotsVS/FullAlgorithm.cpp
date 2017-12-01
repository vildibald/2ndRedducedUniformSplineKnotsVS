#include "stdafx.h"
#include "FullAlgorithm.h"
#include "utils.h"
#include <array>

FullAlgorithm::FullAlgorithm(const InterpolativeMathFunction f)
	: function_(f), xTridiagonals_(),
	  yTridiagonals_(), isParallel_(false)
{
}

Spline FullAlgorithm::Calculate(const KnotVector xVector,
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

void FullAlgorithm::FillDx(Spline& spline)
{
	utils::For(0, static_cast<int>(spline.ColumnsCount()), 1, isParallel_, [&](int j)
	{
		//    for (size_t j = 0; j < spline.ColumnsCount(); ++j) {
		auto& tridiagonal = xTridiagonals_.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto& r2 = tridiagonal.RhsCoeficients()[2];
		const auto& r1 = tridiagonal.RhsCoeficients()[1];
		const auto& r0 = tridiagonal.RhsCoeficients()[0];
		const auto l2 = tridiagonal.Lhs2Coeficients().front();
		const auto l0 = tridiagonal.Lhs0Coeficients().back();

		for (size_t i = 0; i < rightSide.size(); ++i)
		{
			rightSide[i] = r2[i] * spline.Z(i + 2, j) + r1[i] * spline.Z(i + 1, j) + r0[i] * spline.Z(i, j);
		}

		rightSide.front() -= l2 * spline.Dx(0, j);
		rightSide.back() -= l0 * spline.Dx(spline.RowsCount() - 1, j);
		tridiagonal.Solve();

		for (auto i = 0; i < rightSide.size(); ++i)
		{
			spline.SetDx(i + 1, j, rightSide[i]);
		}
	});
}

void FullAlgorithm::FillDy(Spline& spline)
{
	utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i)
	{
		//    for (size_t i = 0; i < spline.RowsCount(); ++i) {
		auto& tridiagonal = yTridiagonals_.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto& r2 = tridiagonal.RhsCoeficients()[2];
		const auto& r1 = tridiagonal.RhsCoeficients()[1];
		const auto& r0 = tridiagonal.RhsCoeficients()[0];
		const auto l2 = tridiagonal.Lhs2Coeficients().front();
		const auto l0 = tridiagonal.Lhs0Coeficients().back();

		for (size_t j = 0; j < rightSide.size(); ++j)
		{
			rightSide[j] = r2[j] * spline.Z(i, j + 2) + r1[j] * spline.Z(i, j + 2) + r0[j] * spline.Z(i, j);
		}

		rightSide.front() -= l2 * spline.Dy(i, 0);
		rightSide.back() -= l0 * spline.Dy(i, spline.ColumnsCount() - 1);
		tridiagonal.Solve();
		for (int j = 0; j < rightSide.size(); ++j)
		{
			spline.SetDy(i, j + 1, rightSide[j]);
		}
	});
}

void FullAlgorithm::FillDxy(Spline& spline)
{
	auto& tridiagonal = xTridiagonals_.Get();
	auto& rightSide = tridiagonal.RightSideBuffer();
	const auto& r2 = tridiagonal.RhsCoeficients()[2];
	const auto& r1 = tridiagonal.RhsCoeficients()[1];
	const auto& r0 = tridiagonal.RhsCoeficients()[0];
	const auto l2 = tridiagonal.Lhs2Coeficients().front();
	const auto l0 = tridiagonal.Lhs0Coeficients().back();
	
	int loop[] = { 0, spline.ColumnsCount() - 1};
	for (auto j : loop)
	{
		for (size_t i = 0; i < rightSide.size(); ++i)
		{
			rightSide[i] = r2[i] * spline.Dy(i + 2, j) + r1[i] * spline.Dy(i + 1, j) + r0[i] * spline.Dy(i, j);
		}
		rightSide.front() -= l2 * spline.Dxy(0, j);
		rightSide.back() -= l0 * spline.Dxy(spline.RowsCount() - 1, j);
		tridiagonal.Solve();
		for (int i = 0; i < rightSide.size(); ++i)
		{
			spline.SetDxy(i + 1, j, rightSide[i]);
		}
	}
}

void FullAlgorithm::FillDyx(Spline& spline)
{
	utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i)
	{
		//    for (size_t i = 0; i < spline.RowsCount(); ++i) {
		auto& tridiagonal = yTridiagonals_.Get();
		auto& rightSide = tridiagonal.RightSideBuffer();
		const auto& r2 = tridiagonal.RhsCoeficients()[2];
		const auto& r1 = tridiagonal.RhsCoeficients()[1];
		const auto& r0 = tridiagonal.RhsCoeficients()[0];
		const auto l2 = tridiagonal.Lhs2Coeficients().front();
		const auto l0 = tridiagonal.Lhs0Coeficients().back();

		for (size_t j = 0; j < rightSide.size(); ++j)
		{
			rightSide[j] = r2[j] * spline.Dx(i, j + 2) + r1[j] * spline.Dx(i, j + 2) + r0[j] * spline.Dx(i, j);
		}

		rightSide.front() -= l2 * spline.Dxy(i, 0);
		rightSide.back() -= l0 * spline.Dxy(i, spline.ColumnsCount() - 1);
		tridiagonal.Solve();
		for (int j = 0; j < rightSide.size(); ++j)
		{
			spline.SetDxy(i, j + 1, rightSide[j]);
		}
	});
}

void FullAlgorithm::Initialize(Spline& spline)
{
	InitializeTridiagonals(spline);
	spline.Initialize(function_);
}

bool FullAlgorithm::IsParallel() const
{
	return isParallel_;
}

void FullAlgorithm::InParallel(const bool value)
{
	isParallel_ = value;
}

void FullAlgorithm::Parallelize(const bool inParallel)
{
	xTridiagonals_.Parallelize(inParallel);
	yTridiagonals_.Parallelize(inParallel);
}

void FullAlgorithm::InitializeTridiagonals(Spline& spline)
{
	xTridiagonals_.GetAll().clear();
	xTridiagonals_.GetAll().clear();
	xTridiagonals_.GetAll().emplace_back(
		Tridiagonal::Factory::CreateFullTridiagonal(spline.X, spline.RowsCount()));
	yTridiagonals_.GetAll().emplace_back(
		Tridiagonal::Factory::CreateFullTridiagonal(spline.Y, spline.ColumnsCount()));
	Parallelize(isParallel_);
}

double FullAlgorithm::ExecutionTime()
{
	return timer_.ExecutionTime();
}

double FullAlgorithm::AllTime()
{
	return timer_.AllTime() +
		xTridiagonals_.GetAll()[0].AllTime() + yTridiagonals_.GetAll()[0].AllTime();
}

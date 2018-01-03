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
	const auto parameterGetter = [&](size_t i, size_t j)
	{
		return spline.Z(i, j);
	};
	const auto derivationGetter = [&](size_t i, size_t j)
	{
		return spline.Dx(i, j);
	};
	auto derivationSetter = [&](size_t i, size_t j, double value)
	{
		spline.SetDx(i, j, value);
	};
	FillD(spline.ColumnsCount(), spline.RowsCount(), xTridiagonals_, parameterGetter,
	      derivationGetter, derivationSetter);
}

void ReducedAlgorithm::FillDy(Spline& spline)
{
	const auto parameterGetter = [&](size_t i, size_t j)
	{
		return spline.Z(j, i);
	};
	const auto derivationGetter = [&](size_t i, size_t j)
	{
		return spline.Dy(j, i);
	};
	auto derivationSetter = [&](size_t i, size_t j, double value)
	{
		spline.SetDy(j, i, value);
	};
	FillD(spline.RowsCount(), spline.ColumnsCount(), yTridiagonals_, parameterGetter,
	      derivationGetter, derivationSetter);
}

void ReducedAlgorithm::FillDxy(Spline& spline)
{
	const auto parameterGetter = [&](size_t i, size_t j)
	{
		return spline.Dy(i, j);
	};
	const auto derivationGetter = [&](size_t i, size_t j)
	{
		return spline.Dxy(i, j);
	};
	auto derivationSetter = [&](size_t i, size_t j, double value)
	{
		spline.SetDxy(i, j, value);
	};

	size_t loop[] = {0, spline.ColumnsCount() - 1};
	for (auto j : loop)
	{
		Solve(spline.ColumnsCount(), spline.RowsCount(), xTridiagonals_,
		      parameterGetter, derivationGetter, derivationSetter, j);
	}
}

void ReducedAlgorithm::FillDyx(Spline& spline)
{
	const auto parameterGetter = [&](size_t i, size_t j)
	{
		return spline.Dx(j, i);
	};
	const auto derivationGetter = [&](size_t i, size_t j)
	{
		return spline.Dxy(j, i);
	};
	auto derivationSetter = [&](size_t i, size_t j, double value)
	{
		spline.SetDxy(j, i, value);
	};
	FillD(spline.RowsCount(), spline.ColumnsCount(), yTridiagonals_, parameterGetter,
	      derivationGetter, derivationSetter);
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
		Tridiagonal::Factory::CreateReducedTridiagonal(spline.X(), spline.RowsCount()));
	yTridiagonals_.GetAll().emplace_back(
		Tridiagonal::Factory::CreateReducedTridiagonal(spline.Y(), spline.ColumnsCount()));

	Parallelize(isParallel_);
}

double ReducedAlgorithm::ExecutionTime()
{
	return timer_.ExecutionTime();
}

double ReducedAlgorithm::AllTime()
{
	return timer_.AllTime() + xTridiagonals_.GetAll()[0].AllTime() + yTridiagonals_.GetAll()[0].
		AllTime();
}

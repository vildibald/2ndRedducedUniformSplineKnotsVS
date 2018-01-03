#include "stdafx.h"
#include "Tridiagonal.h"
#include "utils.h"
#include <algorithm>

std::vector<Tridiagonal>& Tridiagonals::GetAll()
{
	return tridiagonals_;
}

Tridiagonal& Tridiagonals::Get()
{
	if (tridiagonals_.size() == 1)
	{
		return tridiagonals_[0];
	}
	return tridiagonals_[omp_get_thread_num()];
}

void Tridiagonals::Initialize(Tridiagonal tridiagonal)
{
	GetAll().clear();
	GetAll().emplace_back(std::move(tridiagonal));
}

std::vector<double>&
Tridiagonal::ResetBufferAndGet()
{
	auto& buffer = luBuffer_;
	std::fill(buffer.begin(), buffer.end(), 1);
	return buffer;
}

std::vector<double>&
Tridiagonal::Buffer()
{
	return luBuffer_;
}

std::vector<double>&
Tridiagonal::Solve()
{
	auto& buffer = Buffer();
	utils::SolveDeboorTridiagonalSystemBuffered(lhs2Coeficient_,
	                                            lhs1Coeficient_,
	                                            lhs0Coeficient_,
	                                            &rightSideBuffer_.front(), equationCount_,
	                                            &buffer.front());
	return rightSideBuffer_;
}

Tridiagonal::Tridiagonal(const double lhs0Coeficient,
                         const double lhs1Coeficient,
                         const double lhs2Coeficient,
                         const double h,
                         const size_t equationCount,
                         const size_t problemSize)
	: luBuffer_(equationCount),
	  rightSideBuffer_(equationCount),
	  lhs0Coeficient_(lhs0Coeficient),
	  lhs1Coeficient_(lhs1Coeficient),
	  lhs2Coeficient_(lhs2Coeficient),
	  h_(h),
	  equationCount_(equationCount),
	  problemSize_(problemSize)
{
	luBuffer_.assign(equationCount, 0);
	rightSideBuffer_.assign(equationCount, 0);
}


std::vector<double>& Tridiagonal::RightSideBuffer()
{
	return rightSideBuffer_;
}


double Tridiagonal::Lhs0Coeficient() const
{
	return lhs0Coeficient_;
}

double Tridiagonal::Lhs1Coeficient() const
{
	return lhs1Coeficient_;
}

double Tridiagonal::Lhs2Coeficient() const
{
	return lhs2Coeficient_;
}

double Tridiagonal::H() const
{
	return h_;
}

size_t Tridiagonal::EquationCount() const
{
	return equationCount_;
}

size_t Tridiagonal::ProblemSize() const
{
	return problemSize_;
}


Tridiagonal Tridiagonal::Factory::CreateEmptyTridiagonal()
{
	Tridiagonal tridiagonal(0,
	                        0,
	                        0,
	                        0,
	                        0,
	                        0
	);
	return tridiagonal;
}

Tridiagonal
Tridiagonal::Factory::CreateFullTridiagonal(const KnotVector& knotVector,
                                            const size_t numKnots)
{
	const auto numUnknowns = numKnots - 2;
	const auto h = knotVector[1] - knotVector[0];
	Tridiagonal tridiagonal(
		1,
		4,
		1,
		h,
		numUnknowns,
		knotVector.size()
	);
	return tridiagonal;
}

Tridiagonal
Tridiagonal::Factory::CreateReducedTridiagonal(const KnotVector& knotVector,
                                               const size_t numKnots)
{
	const auto remainder = numKnots % 4;
	int numUnknowns = numKnots / 4;
	switch (remainder)
	{
	case 3:
	case 2:
		// Nothing to do.
		break;
	case 1:
	default:
		--numUnknowns;
	}
	const auto h = knotVector[1] - knotVector[0];
	Tridiagonal tridiagonal(
		1,
		-194,
		1,
		h,
		numUnknowns,
		knotVector.size()
	);
	return tridiagonal;
}

double Tridiagonal::AccumulateExecutionTimes(Tridiagonals& tridiagonals)
{
	return 0;
}

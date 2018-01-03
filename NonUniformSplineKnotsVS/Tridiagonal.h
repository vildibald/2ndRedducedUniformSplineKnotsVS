#pragma once

#include "stdafx.h"
#include <vector>
#include "Spline.h"
#include "KnotVector.h"
#include "MultithreadPreparator.h"
#include <numeric>
#include "Timer.h"


class Tridiagonal;

class Tridiagonals;

class Tridiagonals final
{
	std::vector<Tridiagonal> tridiagonals_;

public:
	std::vector<Tridiagonal>& GetAll();

	Tridiagonal& Get();

	void Parallelize(const bool inParallel)
	{
		MultithreadPreparator multithreadPreparator;
		multithreadPreparator.PrepareVector(inParallel, GetAll());
	}

	void Initialize(Tridiagonal tridiagonal);
};

class Tridiagonal final
{
	std::vector<double> luBuffer_{};
	std::vector<double> rightSideBuffer_{};
	double lhs0Coeficient_;
	double lhs1Coeficient_;
	double lhs2Coeficient_;
	double h_;
	size_t equationCount_;
	size_t problemSize_;
	Timer timer_;

	Tridiagonal(double lhs0Coeficient,
	            double lhs1Coeficient,
	            double lhs2Coeficient,
				double h,
	            size_t equationCount,
	            size_t problemSize
	);

public:

	std::vector<double>& Solve();

	std::vector<double>& ResetBufferAndGet();

	std::vector<double>& Buffer();

	double Lhs0Coeficient() const;

	double Lhs1Coeficient() const;

	double Lhs2Coeficient() const;

	double H() const;

	std::vector<double>& RightSideBuffer();

	size_t EquationCount() const;

	size_t ProblemSize() const;

	double ExecutionTime() {
		return timer_.ExecutionTime();
	}

	double AllTime() {
		return timer_.AllTime();
	}

	class Factory final
	{
	public:

		static Tridiagonal
		CreateEmptyTridiagonal();

		static Tridiagonal
		CreateFullTridiagonal(const KnotVector& knotVector, size_t numKnots);

		static Tridiagonal
		CreateReducedTridiagonal(const KnotVector& knotVector, size_t numKnots);

	};

	static double AccumulateAllTimes(Tridiagonals& tridiagonals) {
		return std::accumulate(tridiagonals.GetAll().begin(), tridiagonals.GetAll().end(),
			tridiagonals.GetAll()[0].AllTime(),
			[](double time, Tridiagonal& tridiagonal) {
			return time + tridiagonal.AllTime();
		});
	}

	static double AccumulateExecutionTimes(Tridiagonals& tridiagonals);
};

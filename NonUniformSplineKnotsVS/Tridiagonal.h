#pragma once

#include "stdafx.h"
#include <vector>
#include "Spline.h"
#include "KnotVector.h"


namespace splineknots {
	class Tridiagonal;

	typedef std::vector<Tridiagonal> Tridiagonals;

	class Tridiagonal final {
		std::vector<double> luBuffer;
		std::vector<double> rightSideBuffer;
		std::vector<double> lhs0Coeficients;
		std::vector<double> lhs1Coeficients;
		std::vector<double> lhs2Coeficients;
		std::vector<std::vector<double>> rhsCoeficients;
		size_t numUnknowns;
		size_t problemSize;

		Tridiagonal(std::vector<double> lhs0Coeficients,
			std::vector<double> lhs1Coeficients,
			std::vector<double> lhs2Coeficients,
			std::vector<std::vector<double>> rhsCoeficients,
			size_t numUnknowns,
			size_t problemSize
		);

	public:

		std::vector<double>& Solve();

		std::vector<double>& ResetBufferAndGet();

		std::vector<double>& Buffer();

		const std::vector<double>& GetLhs0Coeficients() const;

		const std::vector<double>& GetLhs1Coeficients() const;

		const std::vector<double>& GetLhs2Coeficients() const;

		const std::vector<std::vector<double>>& GetRhsCoeficients() const;

		std::vector<double>& GetRightSideBuffer();

		size_t GetNumUnknowns() const;

		size_t GetProblemSize() const;

		class Factory final {
		public:
			static Tridiagonal
				CreateFullTridiagonal(const KnotVector& knotVector, size_t numUnknowns);

			static Tridiagonal
				CreateEnhancedTridiagonal(const KnotVector& knotVector, size_t numUnknowns);

		};
	};
}
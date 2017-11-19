#include "stdafx.h"
#include "Tridiagonal.h"
#include "utils.h"
#include <algorithm>


std::vector<double>&
splineknots::Tridiagonal::ResetBufferAndGet() {
	auto& buffer = luBuffer;
	std::fill(buffer.begin(), buffer.end(), 1);
	return buffer;
}

std::vector<double>&
splineknots::Tridiagonal::Buffer() {
	return luBuffer;
}

std::vector<double>&
splineknots::Tridiagonal::Solve() {
	auto& buffer = Buffer();
	utils::SolveTridiagonalSystemBuffered(&lhs0Coeficients.front(),
		&lhs1Coeficients.front(),
		&lhs2Coeficients.front(),
		&rightSideBuffer.front(), numUnknowns,
		&buffer.front());
	return rightSideBuffer;
}

splineknots::Tridiagonal::Tridiagonal(std::vector<double> lhs0Coeficients,
	std::vector<double> lhs1Coeficients,
	std::vector<double> lhs2Coeficients,
	std::vector<std::vector<double>> rhsCoeficients,
	size_t numUnknowns,
	size_t problemSize)
	: lhs0Coeficients(std::move(lhs0Coeficients)),
	lhs1Coeficients(std::move(lhs1Coeficients)),
	lhs2Coeficients(std::move(lhs2Coeficients)),
	rhsCoeficients(std::move(rhsCoeficients)),
	luBuffer(numUnknowns),
	rightSideBuffer(numUnknowns),
	numUnknowns(numUnknowns),
	problemSize(problemSize) {
	luBuffer.assign(numUnknowns, 0);
	rightSideBuffer.assign(numUnknowns, 0);
}


std::vector<double>& splineknots::Tridiagonal::GetRightSideBuffer() {
	return rightSideBuffer;
}


const std::vector<double>& splineknots::Tridiagonal::GetLhs0Coeficients() const {
	return lhs0Coeficients;
}

const std::vector<double>& splineknots::Tridiagonal::GetLhs1Coeficients() const {
	return lhs1Coeficients;
}

const std::vector<double>& splineknots::Tridiagonal::GetLhs2Coeficients() const {
	return lhs2Coeficients;
}

const std::vector<std::vector<double>>& splineknots::Tridiagonal::GetRhsCoeficients() const {
	return rhsCoeficients;
}

size_t splineknots::Tridiagonal::GetNumUnknowns() const {
	return numUnknowns;
}

size_t splineknots::Tridiagonal::GetProblemSize() const {
	return problemSize;
}


splineknots::Tridiagonal
splineknots::Tridiagonal::Factory::CreateFullTridiagonal(const KnotVector& knotVector,
	size_t numKnots) {
	std::vector<double> lhs0Coeficients;
	std::vector<double> lhs1Coeficients;
	std::vector<double> lhs2Coeficients;
	std::vector<double> rhs0Coeficients;
	std::vector<double> rhs1Coeficients;
	std::vector<double> rhs2Coeficients;

	auto numUnknowns = numKnots - 2;
	lhs0Coeficients.reserve(numUnknowns);
	lhs1Coeficients.reserve(numUnknowns);
	lhs2Coeficients.reserve(numUnknowns);
	rhs0Coeficients.reserve(numUnknowns);
	rhs1Coeficients.reserve(numUnknowns);
	rhs2Coeficients.reserve(numUnknowns);

	for (int i = 1; i < numUnknowns + 1; ++i) {
		auto h1 = knotVector[i + 1] - knotVector[i];
		auto h0 = knotVector[i] - knotVector[i - 1];
		lhs2Coeficients.emplace_back(
			h0
		);

		lhs1Coeficients.emplace_back(
			2 * (h0 + h1)
		);

		lhs0Coeficients.emplace_back(
			h1
		);

		rhs2Coeficients.emplace_back(
			3 * h0 / h1
		);

		rhs1Coeficients.emplace_back(
			3 * (h1 - h0)*(h1 + h0) / (h1*h0)
		);

		rhs0Coeficients.emplace_back(
			-3 * h1 / h0
		);
	}

	std::vector<std::vector<double>> rhsCoeficients =
	{
		std::move(rhs0Coeficients),
		std::move(rhs1Coeficients),
		std::move(rhs2Coeficients)
	};
	Tridiagonal tridiagonal(
		std::move(lhs0Coeficients),
		std::move(lhs1Coeficients),
		std::move(lhs2Coeficients),
		std::move(rhsCoeficients),
		numUnknowns,
		knotVector.size()
	);
	return tridiagonal;
}

splineknots::Tridiagonal
splineknots::Tridiagonal::Factory::CreateEnhancedTridiagonal(const KnotVector& knotVector,
	size_t numKnots) {
	std::vector<double> lhs0Coeficients;
	std::vector<double> lhs1Coeficients;
	std::vector<double> lhs2Coeficients;
	std::vector<double> rhs0Coeficients;
	std::vector<double> rhs1Coeficients;
	std::vector<double> rhs2Coeficients;
	std::vector<double> rhs3Coeficients;
	std::vector<double> rhs4Coeficients;

	auto even = numKnots % 2 == 0;
	auto numUnknowns = even ? numKnots / 2 - 2 : numKnots / 2 - 1;

	lhs0Coeficients.reserve(numUnknowns);
	lhs1Coeficients.reserve(numUnknowns);
	lhs2Coeficients.reserve(numUnknowns);
	rhs0Coeficients.reserve(numUnknowns);
	rhs1Coeficients.reserve(numUnknowns);
	rhs2Coeficients.reserve(numUnknowns);
	rhs3Coeficients.reserve(numUnknowns);
	rhs4Coeficients.reserve(numUnknowns);

	auto until = even ? knotVector.size() - 2 : knotVector.size() - 1;
	for (int i = 2; i < until; i += 2) {

		auto h3 = knotVector[i + 2] - knotVector[i + 1];
		auto h2 = knotVector[i + 1] - knotVector[i];
		auto h1 = knotVector[i] - knotVector[i - 1];
		auto h0 = knotVector[i - 1] - knotVector[i - 2];
		auto h11 = h1*h1;
		auto h22 = h2*h2;
		auto h02 = h0*h2;
		auto h13 = h1*h3;
		auto h0123 = h02*h13;
		auto h3Plush2 = h3 + h2;
		auto h2Plush1 = h2 + h1;
		auto h1Plush0 = h1 + h0;

		lhs0Coeficients.emplace_back(
			(h3Plush2) / h0123
		);

		lhs1Coeficients.emplace_back(
			(h1Plush0 / (h02)+h3Plush2 / (h13)*(1 - 4 * h1Plush0*h2Plush1 / (h02))) / (h1*h2)

		);

		lhs2Coeficients.emplace_back(
			h1Plush0 / h0123
		);

		rhs0Coeficients.emplace_back(
			-3 * h3Plush2 / (h0*h0123)
		);

		rhs1Coeficients.emplace_back(
			3 * h3Plush2*h1Plush0 / (h0123*h11)*(2 - (h0 - h1) / h1)
		);

		rhs2Coeficients.emplace_back(
			3 / (h0123)*
			(h3Plush2*h0 / h11 + h1Plush0 / h22*(2 * h3Plush2*(h1 - h2)*h2Plush1 / h22 - h3))
		);

		rhs3Coeficients.emplace_back(
			-3 * h1Plush0*h3Plush2 / (h0123*h22)*(2 + (h2 - h3) / h3)
		);

		rhs4Coeficients.emplace_back(
			3 * h1Plush0 / (h0123*h3)
		);
	}

	std::vector<std::vector<double>> rhsCoeficients =
	{
		std::move(rhs0Coeficients),
		std::move(rhs1Coeficients),
		std::move(rhs2Coeficients),
		std::move(rhs3Coeficients),
		std::move(rhs4Coeficients)
	};
	Tridiagonal tridiagonal(
		std::move(lhs0Coeficients),
		std::move(lhs1Coeficients),
		std::move(lhs2Coeficients),
		std::move(rhsCoeficients),
		numUnknowns,
		knotVector.size()
	);
	return tridiagonal;
}


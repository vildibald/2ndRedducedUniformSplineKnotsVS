#include "stdafx.h"
#include "SplineKnots.h"
#include <iostream>
#include <algorithm>
#include "ComparisonBenchmarkResult.h"
#include <numeric>
#include "StopWatch.h"
#include "MathFunction.h"
#include "FullAlgorithm.h"
#include "ReducedAlgorithm.h"


void EqualityComparison();


KnotVector TestVector(const double from, const double to, const size_t size) {
	KnotVector vector(size);
	const auto nonUniformity = 0.0;
	auto offset1 = (to - from) / (size - 1) + nonUniformity;
	auto offset2 = (to - from) / (size - 1) + nonUniformity;
	const auto offset = (to - from) / (size - 1);
	auto x = from;
	for (size_t i = 0; i < size; ++i) {
		vector[i] = x;
		//        offset = abs(offset - offset1) < nonUniformity/10 ? offset2 : offset1;

		x += offset;
	}
	return vector;
}

ComparisonBenchmarkResult SurfaceBenchmark(const int numIterations, int numKnots,
                                           const bool inParallel = false, bool optimizedLu = true) {
	const MathFunction function = [](double x, double y) {
		return sin(sqrt(x*x + y*y));
	};

	FullAlgorithm full(function);
	ReducedAlgorithm reduced(function);
	full.InParallel(inParallel);
	reduced.InParallel(inParallel);
	numKnots = numKnots % 2 == 0 ? numKnots + 1 : numKnots;

	std::vector<double> calculated_results;
	std::vector<double> full_times;
	full_times.reserve(numIterations);
	std::vector<double> enhanced_times;
	enhanced_times.reserve(numIterations);
	calculated_results.reserve(numIterations * 3);

	const auto vector = TestVector(-2, 2, numKnots);

	for (size_t i = 0; i < numIterations; i++) {
		auto xKnots = vector;
		auto yKnots = vector;
		auto result = full.Calculate(std::move(xKnots), std::move(yKnots));
		//result.Print();
		calculated_results.push_back(result.Dxy(1, 1));
		full_times.push_back(full.ExecutionTime());
	}

	for (size_t i = 0; i < numIterations; i++)
	{
		KnotVector xKnots = vector;
		KnotVector yKnots = vector;
		auto result = reduced.Calculate(std::move(xKnots), std::move(yKnots));
		//result.Print();
		calculated_results.push_back(result.Dxy(1, 1));
		enhanced_times.push_back(reduced.ExecutionTime());
	}

	const auto fullTime = static_cast<double>(std::accumulate(full_times.begin(),
		full_times.end(), 0))
		/ static_cast<double>(numIterations);
	const auto reducedTime = static_cast<double>(std::accumulate(
		enhanced_times.begin(), enhanced_times.end(), 0))
		/ static_cast<double>(numIterations);
	std::cout << "Ignore " << calculated_results[0] << std::endl;
	ComparisonBenchmarkResult benchmarkResult;
	benchmarkResult.Add(fullTime).Add(reducedTime);
	return benchmarkResult;
}

void PrintSurfaceDeboorResult(ComparisonBenchmarkResult& result) {
	std::cout << "Full : " << result[0] << std::endl;
	//    std::cout << "Reduced : " << result.SecondAlg() << std::endl;
	std::cout << "Enhanced : " << result[1] << std::endl;
	std::cout << "Difference F/E: " << result.Ratio(0, 1) << std::endl;
}

void PrintCurveDeboorResult(ComparisonBenchmarkResult& result) {
	std::cout << "Full : " << result[0] << std::endl;
	//    std::cout << "Reduced : " << result.SecondAlg() << std::endl;
	//    std::cout << "Difference F/R: " << result.Ratio() << std::endl;
}

int main() {
	//#ifdef DEBUG
	//    auto result = SurfaceBenchmark(1, 7, false,true);
	//                PrintSurfaceDeboorResult(result);
	//    return 0;
	//#endif
	const auto optimizedTridiagonals = true;
	while (true) {
		//std::cout << clock();
		// Console clear ...
		// ... for Windows,
		system("cls");
		// ... for Linux/Unix.
		//system("clear");
//		std::cout << "1: Instructions benchmark." << std::endl;
		//        std::cout << "2: Spline curve benchmark." << std::endl;
		std::cout << "3: Spline surface benchmark." << std::endl;
		std::cout << "4: Spline surface benchmark (in parallel)." << std::endl;
		std::cout << "5: Equality comparison." << std::endl;
		std::cout << "Q: End program" << std::endl;
		char input;
		std::cin >> input;
		std::cin.get();
		std::cout << std::endl << "---------------" << std::endl;
		unsigned int numIterations;
		unsigned int numKnots;
		ComparisonBenchmarkResult result;
		switch (input) {
//		case '1':
//			std::cout << "Instructions benchmark" << std::endl <<
//				std::endl;
//			MulDivBenchmark();
//			break;
		case '2':
			//                std::cout << "Spline curve benchmark" << std::endl << std::endl;
			//                std::cout << "Enter number of iterations: " << std::endl;
			//                std::cin >> num_iterations;
			//                std::cout << "Enter number of knots: " << std::endl;
			//                std::cin >> num_knots;
			//                std::cin.get();
			//                result = CurveBenchmark(num_iterations, num_knots,
			//                                        optimized_tridiagonals);
			//                PrintCurveDeboorResult(result);
			break;
		case '3':
			std::cout << "Spline surface benchmark" << std::endl << std::endl;
			std::cout << "Enter number of iterations: " << std::endl;
			std::cin >> numIterations;
			std::cout << "Enter number of knots: " << std::endl;
			std::cin >> numKnots;
			std::cin.get();

			result = SurfaceBenchmark(numIterations, numKnots,
				false, optimizedTridiagonals);
			PrintSurfaceDeboorResult(result);
			break;
		case '4':
			std::cout << "Parallel spline surface benchmark" << std::endl <<
				std::endl;
			std::cout << "Enter number of iterations: " << std::endl;
			std::cin >> numIterations;
			std::cout << "Enter number of knots: " << std::endl;
			std::cin >> numKnots;
			std::cin.get();
			result = SurfaceBenchmark(numIterations, numKnots, true, optimizedTridiagonals);
			PrintSurfaceDeboorResult(result);
			break;
		case '5':
			EqualityComparison();
			break;
		case 'q':
		case 'Q':
			return 0;
		default: ;
		}

		std::cout << "===================" << std::endl;
		/*std::cout << "any key: Restart program." << std::endl;
		std::cout << "Q: End program" << std::endl;*/

		system("pause");
	}
}

void EqualityComparison() {
	const MathFunction function = [](double x, double y) {
		return sin(sqrt(x*x + y*y));
	};

	FullAlgorithm full(function);
	ReducedAlgorithm reduced(function);
	full.InParallel(false);
	reduced.InParallel(false);
	const unsigned int numKnots = 9;

	auto vector = TestVector(-20, 20, numKnots);

	auto xKnots = vector;
	auto yKnots = vector;
	auto resultFull = full.Calculate(std::move(xKnots), std::move(yKnots));

	xKnots = vector;
	yKnots = vector;
	auto resultEnhanced = reduced.Calculate(std::move(xKnots), std::move(yKnots));

	auto minDiffDx = (std::numeric_limits<double>::min)(), maxDiffDx = (std::numeric_limits<double>::min)();
	auto minDiffDy = (std::numeric_limits<double>::max)(), maxDiffDy = (std::numeric_limits<double>::min)();
	auto minDiffDxy = (std::numeric_limits<double>::max)(), maxDiffDxy = (std::numeric_limits<double>::min)();
	auto minDiff = (std::numeric_limits<double>::max)(), maxDiff = (std::numeric_limits<double>::min)();

	InterpolativeMathFunction interpolativeMathFunction = function;
	std::cout << "---------- Knot matrix ----------" << std::endl;
	for (unsigned int i = 0; i < resultFull.RowsCount(); ++i) {
		//    for (unsigned int i = numKnots/2; i < numKnots/2+1; ++i) {
		//    for (unsigned int i = 4; i < 5; ++i) {
		std::cout << "Row " << i << " :\n";
		for (unsigned int j = 0; j < resultFull.ColumnsCount(); ++j) {
			minDiffDx = (std::min)(minDiffDx, std::abs(resultFull.Dx(i, j) - resultEnhanced.Dx(i, j)));
			minDiffDy = (std::min)(minDiffDy, std::abs(resultFull.Dy(i, j) - resultEnhanced.Dy(i, j)));
			minDiffDxy = (std::min)(minDiffDxy, std::abs(resultFull.Dxy(i, j) - resultEnhanced.Dxy(i, j)));

			maxDiffDx = (std::max)(maxDiffDx, std::abs(resultFull.Dx(i, j) - resultEnhanced.Dx(i, j)));
			maxDiffDy = (std::max)(maxDiffDy, std::abs(resultFull.Dy(i, j) - resultEnhanced.Dy(i, j)));
			maxDiffDxy = (std::max)(maxDiffDxy, std::abs(resultFull.Dxy(i, j) - resultEnhanced.Dxy(i, j)));

			minDiff = (std::max)({ minDiff, minDiffDx, minDiffDy, minDiffDxy });
			maxDiff = (std::max)({ maxDiff, maxDiffDx, maxDiffDy, maxDiffDxy });

			std::cout << "\nColumn " << j << ":\n"
				<< "M. z: " << interpolativeMathFunction.Z()(vector[i], vector[j]) << '\n'
				<< "M. dx: " << interpolativeMathFunction.Dx()(vector[i], vector[j]) << '\n'
				<< "M. dy: " << interpolativeMathFunction.Dy()(vector[i], vector[j]) << '\n'
				<< "M. dxy: " << interpolativeMathFunction.Dxy()(vector[i], vector[j]) << "\n\n";
			std::cout << "F. z: " << resultFull.Z(i, j) << '\n'
				<< "F. dx: " << resultFull.Dx(i, j) << '\n'
				<< "F. dy: " << resultFull.Dy(i, j) << '\n'
				<< "F. dxy: " << resultFull.Dxy(i, j) << "\n\n";
			std::cout << "R. z: " << resultEnhanced.Z(i, j) << '\n'
				<< "R. dx: " << resultEnhanced.Dx(i, j) << '\n'
				<< "R. dy: " << resultEnhanced.Dy(i, j) << '\n'
				<< "R. dxy: " << resultEnhanced.Dxy(i, j) << '\n';
		}
		std::cout << std::endl;
	}

	std::cout << "-------------------------------" << std::endl;

	std::cout << "Min diff Dx: " << minDiffDx << std::endl;
	std::cout << "Max diff Dx: " << maxDiffDx << std::endl;
	std::cout << "Min diff Dy: " << minDiffDy << std::endl;
	std::cout << "Max diff Dy: " << maxDiffDy << std::endl;
	std::cout << "Min diff Dxy: " << minDiffDxy << std::endl;
	std::cout << "Max diff Dxy: " << maxDiffDxy << std::endl;
	std::cout << "Min diff: " << minDiff << std::endl;
	std::cout << "Max diff: " << maxDiff << std::endl;
}
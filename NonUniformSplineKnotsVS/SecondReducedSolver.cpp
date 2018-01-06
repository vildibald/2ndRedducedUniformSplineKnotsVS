#include "stdafx.h"
#include "SecondReducedSolver.h"


SecondReducedSolver::SecondReducedSolver()
= default;

Tridiagonal SecondReducedSolver::CreateTridiagonal(const KnotVector& knots)
{
	return Tridiagonal::Factory::CreateSecondReducedTridiagonal(knots);
}

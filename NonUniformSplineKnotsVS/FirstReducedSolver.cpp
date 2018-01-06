#include "stdafx.h"
#include "FirstReducedSolver.h"


FirstReducedSolver::FirstReducedSolver()
= default;

Tridiagonal FirstReducedSolver::CreateTridiagonal(const KnotVector& knots)
{
	return Tridiagonal::Factory::CreateFirstReducedTridiagonal(knots);
}

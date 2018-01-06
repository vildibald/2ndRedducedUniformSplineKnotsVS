#include "stdafx.h"
#include "FullSolver.h"


FullSolver::FullSolver()
= default;

Tridiagonal FullSolver::CreateTridiagonal(const KnotVector& knots)
{
	return Tridiagonal::Factory::CreateFullTridiagonal(knots);
}

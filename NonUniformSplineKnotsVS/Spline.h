#pragma once

#include "KnotVector.h"
#include "KnotMatrix.h"

namespace splineknots
{
	class Spline
	{

		KnotMatrix z;
		KnotMatrix dx;
		KnotMatrix dy;
		KnotMatrix dxy;
		KnotVector x;
		KnotVector y;

		Spline();
	public:
		static Spline EmptySpline();

		bool IsEmpty();

		Spline(KnotVector rowVector, KnotVector columnVector);

		size_t RowsCount() const
		{
			return x.size();
		}

		size_t ColumnsCount() const
		{
			return y.size();
		}


		KnotMatrix Z() const
		{
			return z;
		}

		KnotMatrix Dx() const
		{
			return dx;
		}

		KnotMatrix Dy() const
		{
			return dy;
		}

		KnotMatrix Dxy() const
		{
			return dxy;
		}

		const KnotVector& X() const {
			return x;
		}

		const KnotVector& Y() const {
			return y;
		}

		const double X(size_t i) const {
			return x[i];
		}

		const double Y(size_t j) const {
			return y[j];
		}


		double Z(const size_t i, const size_t j) const
		{
			return z[i][j];
		}
		double Dx(const size_t i, const size_t j) const
		{
			return dx[j][i];
		}
		double Dy(const size_t i, const size_t j) const
		{
			return dy[i][j];
		}
		double Dxy(const size_t i, const size_t j) const
		{
			return dxy[i][j];
		}

		void SetZ(const size_t i, const size_t j, const double value)
		{
			z[i][j] = value;
		}
		void SetDx(const size_t i, const size_t j, const double value)
		{
			dx[j][i] = value;
		}
		void SetDy(const size_t i, const size_t j, const double value)
		{
			dy[i][j] = value;
		}
		void SetDxy(const size_t i, const size_t j, const double value)
		{
			dxy[i][j] = value;
		}

		void Print();
	};
}

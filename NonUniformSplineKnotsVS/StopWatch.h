#pragma once
#include <windows.h>
#include <chrono>


//class StopWatch
//{
//	LARGE_INTEGER frequency_;        // ticks per second
//	LARGE_INTEGER t1_, t2_;           // ticks
//public:
//	StopWatch()
//	{
//		// get ticks per second
//		QueryPerformanceFrequency(&frequency_);
//	}
//
//	void Start()
//	{
//		QueryPerformanceCounter(&t1_);
//	}
//
//	void Stop()
//	{
//		QueryPerformanceCounter(&t2_);
//	}
//
//	double EllapsedTime()
//	{
//		return (t2_.QuadPart - t1_.QuadPart) * 1000.0
//			/ frequency_.QuadPart;
//	}
//};

class StopWatch
{
public:
	StopWatch() = default;
	std::chrono::steady_clock::time_point start_;
	std::chrono::steady_clock::time_point stop_;

	void Start()
	{
		 start_ = std::chrono::high_resolution_clock::now();
	}

	void Stop()
	{
		stop_ = std::chrono::high_resolution_clock::now();
	}

	long long EllapsedTime()
	{
		const auto diff = stop_ - start_;
		return std::chrono::duration_cast<std::chrono::microseconds>(diff).count();
	}
};
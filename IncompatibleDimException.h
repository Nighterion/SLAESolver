#pragma once
#include <exception>

#ifndef OOP_LR_2_IncompatibleDimException_h
#define OOP_LR_2_IncompatibleDimException_h

class IncompatibleDimException : public std::exception
{
private:
	int size_A_1;
	int size_A_2;

	int size_B_1;
	int size_B_2;

public:
	IncompatibleDimException(const char* _Exception, int _s_A_1, int _s_A_2, int _s_B_1, int _s_B_2) : exception(_Exception)
	{
		size_A_1 = _s_A_1;
		size_A_2 = _s_A_2;
		size_B_1 = _s_B_1;
		size_B_2 = _s_B_2;
	}
	
	int GetSA1()
	{
		return size_A_1;
	}

	int GetSA2()
	{
		return size_A_2;
	}

	int GetSB1()
	{
		return size_B_1;
	}

	int GetSB2()
	{
		return size_B_2;
	}
};

#endif
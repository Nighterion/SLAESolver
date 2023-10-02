
#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include "BiCGSTAB.h"

void tester(BiCGSTAB& problem)
{
	std::vector<int> threads({ 1,2,3,4,6,8,12,16,32,64,128,256});

	for (int i = 0; i < threads.size(); i++)
	{
		double precision = 0.1;

		omp_set_num_threads(threads.at(i));

		std::cout << std::endl << "threads (" << threads.at(i) << ")" << std::endl;

		double timeStart = omp_get_wtime();
		try
		{
			problem.solve(precision);
		}
		catch (std::exception exception)
		{
			std::cout << "Exception caught :" << std::endl << "''" << exception.what() << "''";
			return;
		}
		double timeEnd = omp_get_wtime();

		std::cout << std::endl << "------------ Summary ------------" << std::endl << "Precision : " << problem.evaluateSolution() << std::endl;

		std::cout << std::endl << "Time elapsed : " << timeEnd - timeStart << std::endl;

		std::cout << std::endl << "Iterations elapsed : " << problem.getIterCount() << std::endl << "---------------------------------" << std::endl;
	}
}

void printvector(std::vector<double> p_vec)
{
	std::cout << "-*-*-*-*-*-*-*-" << std::endl;

    for (int i = 0; i < p_vec.size(); i++)
    {
        std::cout << std::setw(4) << i + 1 << " ->    " << p_vec.at(i) << std::endl;
    }

	std::cout << "-*-*-*-*-*-*-*-" << std::endl << std::endl;
}

int main()
{
	std::string matrixFName;

	std::string vectorFName;

	srand(time(NULL));

    omp_set_num_threads(4);

	std::cout << "Input matrix file name >>>  ";
	std::cin >> matrixFName;
	std::cout << std::endl;

	std::cout << "Input vector file name >>>  ";
	std::cin >> vectorFName;
	std::cout << std::endl;

	BiCGSTAB problem;

	try
	{
		problem.setA(matrixFName);
		problem.setb(vectorFName);
	}
	catch (std::exception exception)
	{
		std::cout << "Exception caught :" << std::endl << "''" << exception.what() << "''";
		return 1;
	}

	std::cout << "Current problem dimensions :  " << problem.getSize() << std::endl << std::endl;

	//tester(problem);

	//return 0;

	while (threads != 0)
	{
		std::string option;

		double precision;

		std::cout << std::endl << "Input thread number >>>  ";


		std::cin >> threads;

		std::cout << std::endl << "Input precision >>>  ";

		std::cin >> precision;



		omp_set_num_threads(threads);

        

		double timeStart = omp_get_wtime();
		try
		{
			problem.solve(precision);
		}
		catch (std::exception exception)
		{
			std::cout << "Exception caught :" << std::endl << "''" << exception.what() << "''";
			return 1;
		}
		double timeEnd = omp_get_wtime();

        std:: cout << std::endl << "------------ Summary ------------" << std::endl << "Precision : " << problem.evaluateSolution() << std::endl;

        std::cout << std::endl << "Time elapsed : " << timeEnd - timeStart << std::endl;

        std::cout << std::endl << "Iterations elapsed : " << problem.getIterCount() << std::endl << "---------------------------------" << std::endl;

		bool badInput = true;

		while (badInput)
		{
			std::cout << std::endl << "Show answer?(Y/N) >>>  ";
			std::cin >> option;

			if (option == "Y" || option == "y")
			{
				auto res = problem.getSolution();
				printvector(res);
				badInput = false;
			}
			else if (option == "N" || option == "n")
			{
				badInput = false;
			}
		}

		badInput = true;

		while (badInput)
		{
			std::cout << std::endl << "Exit?(Y/N) >>>  ";
			std::cin >> option;

			if (option == "Y" || option == "y")
			{
				return 0;
			}
			else if (option == "N" || option == "n")
			{
				badInput = false;
			}
		}
	}
	return 0;

}
#include "ff.h"

#include <cmath>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

int main()
{
	float aSide = 1.0f;
	float bSide = 1.0f;
	float cSide = 1.0f;

	float poly0[4][3] =
	{
		{ 0.0f, 0.0f, 0.0f },
		{ -aSide, 0.0f, 0.0f },
		{ -aSide, 0.0f, cSide },
		{ 0.0f, 0.0f, cSide }
	};

	float poly1[4][3] =
	{
		{ 0.0f, 0.0f, 0.0f },
		{ 0.0f, 0.0f, cSide },
		{ bSide, 0.0f, cSide },
		{ bSide, 0.0f, 0.0f }
	};

	unsigned int angleStep = 512;
	unsigned int sideStep = 5;

	std::ofstream ffFile;
	ffFile.open("caltechFF.txt", std::ios_base::out);
	std::stringstream ssFF;

	ssFF << std::fixed;
	for (unsigned int i = 0; i <= angleStep; ++i)
	{
		for (unsigned int j = 1; j <= sideStep; ++j)
		{
			float aSideVar = aSide * static_cast<float>(j) / sideStep;
			poly0[1][0] = poly0[2][0] = static_cast<float>(aSideVar * cos(static_cast<float>(i)* M_PI / static_cast<float>(angleStep)));
			poly0[1][1] = poly0[2][1] = static_cast<float>(aSideVar * sin(static_cast<float>(i)* M_PI / static_cast<float>(angleStep)));
			ssFF << std::setprecision(8) << FormFactorf(poly1, 4, poly0, 4) << "\t";
		}
		ssFF << "\n";
	}

	ffFile << ssFF.str();
	ffFile.close();

	return 0;
}

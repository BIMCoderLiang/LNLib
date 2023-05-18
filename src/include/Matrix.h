#pragma once

#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib
{
	class LNLIB_EXPORT Matrix
	{
	public:
	
		template<typename T>
		static std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>>& matrix)
		{
			std::vector<std::vector<T>> result;
			std::vector<T> temp;

			for (int i = 0; i < static_cast<int>(matrix[0].size()); i++)
			{
				for (int j = 0; j < static_cast<int>(matrix.size()); j++)
				{
					temp.emplace_back(matrix[j][i]);
				}
				result.emplace_back(temp);
				temp.erase(temp.begin(), temp.end());
			}

			return result;
		}
	};
}


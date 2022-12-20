#include <iostream>
#include <vector>
#include <algorithm>

std::vector<std::vector<uint32_t>> getSubsets(std::vector<uint32_t> &nums)
{

    std::vector<std::vector<uint32_t>> subset;
    std::vector<uint32_t> empty;
    subset.push_back(empty);

    for (int i = 0; i < nums.size(); i++)
    {
        std::vector<std::vector<uint32_t>> subsetTemp = subset;

        for (int j = 0; j < subsetTemp.size(); j++)
        {
            subsetTemp[j].push_back(nums[i]);
        }

        for (int j = 0; j < subsetTemp.size(); j++)
        {
            subset.push_back(subsetTemp[j]);
        }
    }
    return subset;
}

int main()
{
    std::vector<uint32_t> nums = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    std::vector<std::vector<uint32_t>> subset = getSubsets(nums);

    for (int i = 0; i < subset.size(); i++)
    {
        for (int j = 0; j < subset[i].size(); j++)
        {
            std::cout << subset[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
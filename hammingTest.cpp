#include <cstdint>
#include <iostream>
#include <vector>
#include <chrono>

uint32_t hammingDistance(const uint32_t bit_a, const uint32_t bit_b)
{
    return __builtin_popcount(bit_a ^ bit_b);
}

std::string binaryRep(const uint32_t bit, const uint32_t bitlength)
{
    std::string s;
    s.reserve(bitlength);
    s = "";

    for (uint32_t i = 0; i < bitlength; i++)
    {
        s = ((bit & (1 << i)) ? "1" : "0") + s;
    }
    return s;
}

std::vector<uint32_t> generateOne(uint32_t length, uint32_t d, uint32_t n)
{
    std::vector<uint32_t> adjacencies;
    uint32_t mask = (1u << d) - 1;

    // Find all numbers with Hamming distance equal to or greater than d from n
    for (uint32_t i = d; i < length; i++)
    {
        while (mask < (1u << length))
        {

            adjacencies.push_back(mask ^ n);

            // Shift the mask to the next permutation
            uint32_t t = mask | (mask - 1);
            mask = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(mask) + 1));
        }
    }
    return adjacencies;
}

std::vector<std::vector<uint32_t>> generateAll_cheap(uint32_t length, uint32_t d)
{
    std::vector<std::vector<uint32_t>> adjacencies;
    adjacencies.reserve(1u << length);
    for (uint32_t i = 0; i < (1u << length); i++)
    {
        adjacencies.push_back(generateOne(length, d, i));
    }
    return adjacencies;
}


std::vector<std::vector<uint32_t>> generateAllFast(uint32_t length, uint32_t d)
{
    const uint32_t max = 1u << length;
    std::vector<std::vector<uint32_t>> adjacencies;
    adjacencies.reserve(max);
    for (uint32_t i = 0; i < max; i++)
    {
        adjacencies.push_back(std::vector<uint32_t>());
        for (uint32_t j = d; j <= length; j++)
        {
            uint32_t mask = (1u << j) - 1;
            while (mask < max)
            {

                adjacencies[i].push_back(mask ^ i);

                // Shift the mask to the next permutation
                uint32_t t = mask | (mask - 1);
                mask = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(mask) + 1));
            }
        }
    }
    return adjacencies;
}

std::vector<std::vector<uint32_t>> generateAllSlow(uint32_t length, uint32_t d)
{
    std::vector<std::vector<uint32_t>> adjacencies;
    adjacencies.reserve(1u << length);
    for (uint32_t i = 0; i < adjacencies.size(); i++)
    {
        adjacencies.push_back(std::vector<uint32_t>());
        for (uint32_t j = i + 1; j < adjacencies.size(); j++)
        {
            if (hammingDistance(i, j) >= d)
            {
                adjacencies[i].push_back(j);
                adjacencies[j].push_back(i);
            }
        }
    }
    return adjacencies;
}

int main()
{

    uint32_t length = 12;
    uint32_t d = 2;

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<uint32_t>> adjacencies = generateAllFast(length, d);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Fast: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();

    adjacencies = generateAllSlow(length, d);

    end = std::chrono::high_resolution_clock::now();

    std::cout << "Slow: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}
#include <cstdint>
#include <vector>
#include <iostream>

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

uint32_t hammingDistance(const uint32_t bit_a, const uint32_t bit_b)
{
    return __builtin_popcount(bit_a ^ bit_b);
}

std::vector<uint32_t> findAdjacency(uint32_t d, uint32_t number, uint32_t bitlength)
{
    std::vector<uint32_t> adjacency;
    uint32_t m = (1 << d) - 1; // create bit mask with d 1's
    for (uint32_t c = 0; c < (1 << d); c++)
    {                                                         // iterate over all combinations of d bits
        uint32_t flipped = number ^ (c & m);                  // flip d bits specified by c
        uint32_t distance = hammingDistance(number, flipped); // calculate distance using hammingDistance function
        if (distance >= d)
        { // check if distance is at least d
            adjacency.push_back(flipped);
        }
    }
    return adjacency;
}

int main()
{
    uint32_t d = 4;
    uint32_t bitlength = 5;
    uint32_t number = 0;
    std::vector<uint32_t> adjacency = findAdjacency(d, number, bitlength);

    for (uint32_t i = 0; i < adjacency.size(); i++)
    {
        std::cout << binaryRep(adjacency[i], bitlength) << std::endl;
    }

    return 0;
}
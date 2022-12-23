#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

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
// std::bitset::count might be better

class graph
{
public:
    std::vector<uint32_t> nodes;
    uint32_t bitlength;
    uint32_t distance;

    // Default constructor
    graph()
    {
    }

    // Constructor
    graph(uint32_t bitlength, uint32_t distance):bitlength(bitlength), distance(distance)
    {
        for (uint32_t i = 0; i < (1 << bitlength); i++)
        {
            nodes.push_back(i);
        }
    }

    // Copy constructor
    graph(const graph &g)
    {
        this->bitlength = g.bitlength;
        this->distance = g.distance;
        this->nodes = g.nodes;
    }

    // Destructor
    ~graph()
    {
        nodes.clear();
    }

    static bool isSameNodes(const std::vector<uint32_t> &nodes_a, const std::vector<uint32_t> &nodes_b)
    {

        std::vector<uint32_t> nodes_b_copy = nodes_b;

        if (nodes_a.size() != nodes_b_copy.size())
        {
            return false;
        }

        for (uint32_t i = 0; i < nodes_a.size(); i++)
        {
            bool found = false;
            for (uint32_t j = 0; j < nodes_b_copy.size(); j++)
            {
                if (nodes_a[i] == nodes_b_copy[j])
                {
                    nodes_b_copy.erase(nodes_b_copy.begin() + j);
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                return false;
            }
        }

        return true;
    }

    static void removeDuplicateNodes(std::vector<std::vector<uint32_t>> &nodes)
    {
        for (uint32_t i = 0; i < nodes.size(); i++)
        {
            for (uint32_t j = i + 1; j < nodes.size(); j++)
            {
                if (isSameNodes(nodes[i], nodes[j]))
                {
                    nodes.erase(nodes.begin() + j);
                    j--;
                }
            }
        }
    }

    std::vector<std::vector<uint32_t>> getSubsets(const std::vector<uint32_t> &nums)
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

    bool checkClique(const std::vector<uint32_t> &clique)
    {
        for (uint32_t i = 0; i < clique.size(); i++)
        {
            for (uint32_t j = i + 1; j < clique.size(); j++)
            {
                if (hammingDistance(clique[i], clique[j]) < distance)
                {
                    return false;
                }
            }
        }
        return true;
    }

    bool checkMaxmialClique(std::vector<uint32_t> &clique)
    {
        for (uint32_t i = 0; i < nodes.size(); i++)
        {
            if (hammingDistance(clique[0], i) >= distance && std::find(clique.begin(), clique.end(), nodes[i]) == clique.end())
            {
                clique.push_back(nodes[i]);
                if (checkClique(clique))
                {
                    clique.pop_back();
                    return false;
                }
                clique.pop_back();
            }
        }
        return true;
    }

    uint32_t findMaximalCliqueBruteForce()
    {
        std::vector<std::vector<uint32_t>> subsets = getSubsets(nodes);

        for (uint32_t i = 0; i < subsets.size(); i++)
        {
            if (!checkClique(subsets[i]) && !checkMaxmialClique(subsets[i]))
            {
                subsets.erase(subsets.begin() + i);
                i--;
            }
        }

        uint32_t max = 0;
        for (uint32_t i = 0; i < subsets.size(); i++)
        {
            if (subsets[i].size() > max)
            {
                max = subsets[i].size();
            }
        }

        return max;
    }

    uint32_t findMaximalCliqueHeuristicBFS()
    {
        std::vector<std::vector<uint32_t>> cliques;
        std::vector<uint32_t> visited(nodes.size(), 0);
        std::vector<uint32_t> current_clique;
        std::vector<uint32_t> queue = {0};
        visited[0] = 1;
        current_clique.push_back(nodes[0]);
        
        while (!queue.empty())
        {
            uint32_t current = queue.back();
            queue.pop_back();
            for (int i = 0; i < nodes.size(); i++)
            {
                if (visited[i] == 0 && hammingDistance(current, nodes[i]) >= distance)
                {
                    bool is_connected_to_all = true;
                    for (const uint32_t &node : current_clique)
                    {
                        if (hammingDistance(nodes[i], node) < distance)
                        {
                            is_connected_to_all = false;
                            break;
                        }
                    }
                    if (is_connected_to_all)
                    {
                        visited[i] = 1;
                        current_clique.push_back(nodes[i]);
                        queue.push_back(nodes[i]);
                    }
                }
            }
            if (queue.empty())
            {
                cliques.push_back(current_clique);
                current_clique.clear();
                for (int i = 0; i < nodes.size(); i++)
                {
                    if (visited[i] == 0)
                    {
                        queue.push_back(i);
                        visited[i] = 1;
                        current_clique.push_back(nodes[i]);
                        break;
                    }
                }
            }
        }

        uint32_t max = 0;

        for (int i = 0; i < cliques.size(); i++)
        {
            if (cliques[i].size() > max)
            {
                max = cliques[i].size();
            }
        }

        return max;
    }

    uint32_t findMaximalCliqueBronKerboschSimple()
    {
        std::vector<uint32_t> R;
        std::vector<uint32_t> P = nodes;
        std::vector<uint32_t> X;
        std::vector<std::vector<uint32_t>> maximalCliques;

        findMaximalCliqueBronKerboschSimple(R, P, X, maximalCliques);

        uint32_t max = 0;
        for (uint32_t i = 0; i < maximalCliques.size(); i++)
        {
            if (maximalCliques[i].size() > max)
            {
                max = maximalCliques[i].size();
            }
        }

        return max;
    }

    void findMaximalCliqueBronKerboschSimple(std::vector<uint32_t> &R, std::vector<uint32_t> &P, std::vector<uint32_t> &X, std::vector<std::vector<uint32_t>> &maximalCliques)
    {
        if (P.empty() && X.empty())
        {
            maximalCliques.push_back(R);
        }
        else
        {
            for (uint32_t i = 0; i < P.size(); i++)
            {
                std::vector<uint32_t> newR = R;
                newR.push_back(P[i]);
                std::vector<uint32_t> newP;
                std::vector<uint32_t> newX;
                for (uint32_t j = 0; j < P.size(); j++)
                {
                    if (hammingDistance(P[i], P[j]) >= distance)
                    {
                        newP.push_back(P[j]);
                    }
                }
                for (uint32_t j = 0; j < X.size(); j++)
                {
                    if (hammingDistance(P[i], X[j]) >= distance)
                    {
                        newX.push_back(X[j]);
                    }
                }
                findMaximalCliqueBronKerboschSimple(newR, newP, newX, maximalCliques);
                for (uint32_t j = 0; j < P.size(); j++)
                {
                    if (P[i] == P[j])
                    {
                        P.erase(P.begin() + j);
                        break;
                    }
                }
                X.push_back(P[i]);
            }
        }
    }

    uint32_t findMaximalCliqueBronKerboschPivot()
    {
        std::vector<uint32_t> R;
        std::vector<uint32_t> P = nodes;
        std::vector<uint32_t> X;
        std::vector<std::vector<uint32_t>> maximalCliques;

        findMaximalCliqueBronKerboschPivot(R, P, X, maximalCliques);

        uint32_t max = 0;
        for (uint32_t i = 0; i < maximalCliques.size(); i++)
        {
            if (maximalCliques[i].size() > max)
            {
                max = maximalCliques[i].size();
            }
        }

        return max;
    }

    uint32_t findMaximalCliqueBronKerboschPivot(std::vector<uint32_t> &R, std::vector<uint32_t> &P, std::vector<uint32_t> &X, std::vector<std::vector<uint32_t>> &maximalCliques)
    {
        // TODO
        return 0;
    }
};

void findMaximalClique(const uint32_t n, const uint32_t d, const uint32_t m)
{

    using clock_t = std::chrono::high_resolution_clock;
    using timepoint_t = std::chrono::time_point<clock_t>;

    timepoint_t start = clock_t::now();

    graph g(n, d);

    uint32_t max = 0;

    switch (m)
    {
    case 1:
        max = g.findMaximalCliqueBruteForce();
        break;
    case 2:
        max = g.findMaximalCliqueHeuristicBFS();
        break;
    case 3:
        max = g.findMaximalCliqueBronKerboschSimple();
        break;
    case 4:
        max = g.findMaximalCliqueBronKerboschPivot();
        break;
    }

    timepoint_t end = clock_t::now();

    double elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0;

    std::cout << "Maximum Clique Size: " << max << std::endl;
    std::cout << "Elapsed Time: " << elapsed_seconds << "s" << std::endl;
}

int main()
{

    std::cout << "Enter n and d to find the maximal clique of a graph with n nodes and d distance." << std::endl;
    std::cout << "Enter m for the method to use." << std::endl;
    std::cout << "Enter 1 for brute force." << std::endl;
    std::cout << "Enter 2 for Heuristic Breadth-First Search." << std::endl;
    std::cout << "Enter 3 for Simple Bron-Kerbosch." << std::endl;
    std::cout << "Enter 4 for Pivot Bron-Kerbosch." << std::endl;
    std::cout << "Enter unvalid numbers to exit." << std::endl;

    while (true)
    {

        std::string temp;

        std::cout << "Enter n: ";
        std::getline(std::cin, temp);

        if (stoi(temp) < 1)
        {
            break;
        }
        uint32_t n = stoi(temp);

        std::cout << "Enter d: ";
        std::getline(std::cin, temp);

        if (stoi(temp) < 1)
        {
            break;
        }
        uint32_t d = stoi(temp);

        std::cout << "Enter m: ";
        std::getline(std::cin, temp);

        if (stoi(temp) < 1 || stoi(temp) > 4)
        {
            break;
        }
        uint32_t m = stoi(temp);

        findMaximalClique(n, d, m);
    }

    return 0;
}
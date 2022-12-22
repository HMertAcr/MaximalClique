#include <iostream>
#include <string>
#include <vector>
#include <queue>
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

class graph
{
public:
    std::vector<uint32_t> nodes;
    uint32_t bitlength;

    // Default constructor
    graph()
    {
    }

    // Constructor
    graph(uint32_t bitlength)
    {
        this->bitlength = bitlength;

        for (uint32_t i = 0; i < (1 << bitlength); i++)
        {
            nodes.push_back(i);
        }
    }

    // Copy constructor
    graph(const graph &g)
    {
        this->bitlength = g.bitlength;
        this->nodes = g.nodes;
    }

    // Destructor
    ~graph()
    {
        nodes.clear();
    }

    static bool isSameNodes(const std::vector<uint32_t> &nodes_a, const std::vector<uint32_t> &nodes_b)
    {

        auto tempNodes_b = nodes_b;

        if (nodes_a.size() != tempNodes_b.size())
        {
            return false;
        }

        for (uint32_t i = 0; i < nodes_a.size(); i++)
        {
            bool found = false;
            for (uint32_t j = 0; j < tempNodes_b.size(); j++)
            {
                if (nodes_a[i] == tempNodes_b[j])
                {
                    tempNodes_b.erase(tempNodes_b.begin() + j);
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

    // TODO: Optimize with adjacency list
    bool checkClique(std::vector<uint32_t> &clique, const uint32_t d)
    {
        for (uint32_t i = 0; i < clique.size(); i++)
        {
            for (uint32_t j = i + 1; j < clique.size(); j++)
            {
                if (hammingDistance(clique[i], clique[j]) < d)
                {
                    return false;
                }
            }
        }
        return true;
    }

    // TODO: Optimize with adjacency list
    bool checkMaxmialClique(std::vector<uint32_t> &clique, const uint32_t d)
    {
        for (uint32_t i = 0; i < nodes.size(); i++)
        {
            if (std::find(clique.begin(), clique.end(), nodes[i]) == clique.end())
            {
                clique.push_back(nodes[i]);
                if (checkClique(clique, d))
                {
                    return false;
                }
                clique.pop_back();
            }
        }
        return true;
    }

    uint32_t findMaximalCliqueBruteForce(const uint32_t d)
    {
        std::vector<std::vector<uint32_t>> subsets = getSubsets(nodes);

        for (uint32_t i = 0; i < subsets.size(); i++)
        {
            if (!checkClique(subsets[i], d) && !checkMaxmialClique(subsets[i], d))
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
    
    //deleted
    uint32_t findMaximalCliqueApprox(const uint32_t d)
    {
        std::vector<std::vector<uint32_t>> clique;

        for (uint32_t i = 0; i < nodes.size(); i++)
        {
            clique.push_back(std::vector<uint32_t>());
            clique[i].push_back(nodes[i]);
            for (uint32_t j = 0; j < nodes.size(); j++)
            {
                bool canAdd = true;
                for (uint32_t k = 0; k < clique[i].size(); k++)
                {
                    if (hammingDistance(nodes[j], clique[i][k]) < d)
                    {
                        canAdd = false;
                        break;
                    }
                }
                if (canAdd)
                {
                    clique[i].push_back(nodes[j]);
                }
            }
        }

        //not needed
        graph::removeDuplicateNodes(clique);

        uint32_t max = 0;
        for (uint32_t i = 0; i < clique.size(); i++)
        {
            if (clique[i].size() > max)
            {
                max = clique[i].size();
            }
        }

        return max;
    }

    // TODO speed up by implementing an adjacency list that calculates adjacencies with dynamic programming and not brute force
    uint32_t findMaximalCliqueHeuristicBFS(const uint32_t d)
    {
        std::vector<std::vector<uint32_t>> cliques;
        std::vector<uint32_t> visited(nodes.size(), 0);
        std::vector<uint32_t> current_clique;
        uint32_t start = 0;
        std::vector<uint32_t> queue = {start};
        visited[start] = 1;
        current_clique.push_back(start);
        while (!queue.empty())
        {
            uint32_t current = queue.back();
            queue.pop_back();
            for (int i = 0; i < nodes.size(); i++)
            {
                if (visited[i] == 0 && hammingDistance(current, nodes[i]) >= d)
                {
                    bool is_connected_to_all = true;
                    for (const uint32_t &node : current_clique)
                    {
                        if (hammingDistance(nodes[i], node) < d)
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

    uint32_t findMaximalCliqueBronKerboschSimple(const uint32_t d)
    {
        std::vector<uint32_t> R;
        std::vector<uint32_t> P = nodes;
        std::vector<uint32_t> X;
        std::vector<std::vector<uint32_t>> maximalCliques;

        findMaximalCliqueBronKerboschSimple(R, P, X, d, maximalCliques);

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

    void findMaximalCliqueBronKerboschSimple(std::vector<uint32_t> R, std::vector<uint32_t> P, std::vector<uint32_t> X, const uint32_t d, std::vector<std::vector<uint32_t>> &maximalCliques)
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
                    if (hammingDistance(P[i], P[j]) >= d)
                    {
                        newP.push_back(P[j]);
                    }
                }
                for (uint32_t j = 0; j < X.size(); j++)
                {
                    if (hammingDistance(P[i], X[j]) >= d)
                    {
                        newX.push_back(X[j]);
                    }
                }
                findMaximalCliqueBronKerboschSimple(newR, newP, newX, d, maximalCliques);
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

    uint32_t findMaximalCliqueBronKerboschPivot(const uint32_t d)
    {
        std::vector<uint32_t> R;
        std::vector<uint32_t> P = nodes;
        std::vector<uint32_t> X;
        std::vector<std::vector<uint32_t>> maximalCliques;

        findMaximalCliqueBronKerboschPivot(R, P, X, d, maximalCliques);

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

    uint32_t findMaximalCliqueBronKerboschPivot(std::vector<uint32_t> R, std::vector<uint32_t> P, std::vector<uint32_t> X, const uint32_t d, std::vector<std::vector<uint32_t>> &maximalCliques)
    {
        // TODO
        return 0;
    }
};

void findMaximalClique(const uint32_t n, const uint32_t d, const uint32_t m)
{

    auto start = std::chrono::high_resolution_clock::now();

    graph g(n);

    uint32_t max = 0;

    switch (m)
    {
    case 1:
        max = g.findMaximalCliqueBruteForce(d);
        break;
    case 2:
        max = g.findMaximalCliqueHeuristicBFS(d);
        break;
    case 3:
        max = g.findMaximalCliqueBronKerboschSimple(d);
        break;
    case 4:
        max = g.findMaximalCliqueBronKerboschPivot(d);
        break;
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Maximum Clique Size: " << max << std::endl;
    std::cout << "Elapsed Time: " << elapsed.count() << "s" << std::endl;
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
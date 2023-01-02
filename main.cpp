#include <iostream>
#include <string>
#include <sstream>
#include <bitset>
#include <vector>
#include <unordered_set>
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
    return std::bitset<32>(bit_a ^ bit_b).count();
}

template <typename Container>
uint32_t getLargestListSize(const Container &cliques)
{
    return std::max_element(cliques.begin(), cliques.end(), [](const auto &a, const auto &b)
                            { return a.size() < b.size(); })
        ->size();
}

class Adjacency
{
public:
    uint32_t bitlength;
    uint32_t distance;
    bool isAdjacencyCreated = false;

    // Default constructor
    Adjacency()
    {
    }

    // Constructor
    Adjacency(uint32_t bitlength, uint32_t distance) : distance(distance), bitlength(bitlength)
    {
        uint32_t s = (1 << bitlength);
        for (uint32_t i = 0; i < s; i++)
        {
            adjacencies.push_back(std::vector<uint32_t>());
            adjacencies[i].reserve(s);
        }
    }

    // Copy constructor
    Adjacency(const Adjacency &a)
    {
        this->distance = a.distance;
        this->adjacencies = a.adjacencies;
    }

    // Destructor
    ~Adjacency()
    {
        adjacencies.clear();
    }

    void createAdjacencies()
    {
        if (!isAdjacencyCreated)
        {
            for (uint32_t i = 0; i < adjacencies.size(); i++)
            {
                for (uint32_t j = i + 1; j < adjacencies.size(); j++)
                {
                    if (hammingDistance(i, j) >= distance)
                    {
                        adjacencies[i].push_back(j);
                        adjacencies[j].push_back(i);
                    }
                }
            }
            isAdjacencyCreated = true;
        }
    }

    // Returns the set of d-adjacent nodes for the given node 'x'.
    const std::vector<uint32_t> &getAdjacencies(uint32_t x)
    {
        if (!isAdjacencyCreated)
        {
            if (adjacencies[x].size() == 0)
            {
                for (uint32_t i = 0; i < adjacencies.size(); i++)
                {
                    if (hammingDistance(x, i) >= distance)
                    {
                        adjacencies[x].push_back(i);
                    }
                }
            }
        }
        return adjacencies[x];
    }

    // Returns true if nodes 'a' and 'b' are d-adjacent, false otherwise.
    bool areAdjacent(uint32_t a, uint32_t b) const
    {
        return hammingDistance(a, b) >= distance;
    }

private:
    std::vector<std::vector<uint32_t>> adjacencies;
};

class graph
{
public:
    std::vector<uint32_t> nodes;
    Adjacency adjacency;
    uint32_t bitlength;
    uint32_t distance;

    // Default constructor
    graph()
    {
    }

    // Constructor
    graph(uint32_t bitlength, uint32_t distance) : bitlength(bitlength), distance(distance)
    {
        for (uint32_t i = 0; i < (1u << bitlength); i++)
        {
            nodes.push_back(i);
        }

        adjacency = Adjacency(bitlength, distance);
    }

    // Copy constructor
    graph(const graph &g)
    {
        this->bitlength = g.bitlength;
        this->distance = g.distance;
        this->nodes = g.nodes;
        this->adjacency = g.adjacency;
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
                if (adjacency.areAdjacent(clique[i], clique[j]))
                {
                    return false;
                }
            }
        }
        return true;
    }

    bool checkMaxmialClique(std::vector<uint32_t> &clique)
    {
        const std::vector<uint32_t> &adjacencies = adjacency.getAdjacencies(clique[0]);
        for (const uint32_t &node : adjacencies)
        {
            if (std::find(clique.begin(), clique.end(), node) == clique.end())
            {
                clique.push_back(node);
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

        adjacency.createAdjacencies();

        std::vector<std::vector<uint32_t>> subsets = getSubsets(nodes);

        for (uint32_t i = 0; i < subsets.size(); i++)
        {
            if (!checkClique(subsets[i]) && !checkMaxmialClique(subsets[i]))
            {
                subsets.erase(subsets.begin() + i);
                i--;
            }
        }

        return getLargestListSize(subsets);
    }

    uint32_t findMaximalCliqueHeuristicBFS()
    {
        std::vector<std::unordered_set<uint32_t>> cliques;

        for (const uint32_t &start_node : nodes)
        {
            std::queue<uint32_t> node_queue({start_node});
            std::unordered_set<uint32_t> clique_nodes({start_node});

            while (!node_queue.empty())
            {
                uint32_t current_node = node_queue.front();
                node_queue.pop();

                for (const uint32_t &neighbor : adjacency.getAdjacencies(current_node))
                {
                    if (clique_nodes.count(neighbor))
                    {
                        continue;
                    }
                    bool is_adjacent_to_all = true;
                    for (const uint32_t &clique_node : clique_nodes)
                    {
                        if (!adjacency.areAdjacent(neighbor, clique_node))
                        {
                            is_adjacent_to_all = false;
                            break;
                        }
                    }
                    if (is_adjacent_to_all)
                    {
                        clique_nodes.insert(neighbor);
                        node_queue.push(neighbor);
                    }
                }
            }

            cliques.push_back(clique_nodes);
        }

        return getLargestListSize(cliques);
    }

    uint32_t findMaximalCliqueGreedy()
    {
        std::vector<std::vector<uint32_t>> cliques;
        std::vector<uint32_t> remaining_nodes = nodes;

        while (!remaining_nodes.empty())
        {
            uint32_t max_degree_vertex = remaining_nodes[0];

            uint32_t max_degree = adjacency.getAdjacencies(max_degree_vertex).size();

            for (const uint32_t &v : remaining_nodes)
            {
                uint32_t degree = adjacency.getAdjacencies(v).size();
                if (degree > max_degree)
                {
                    max_degree_vertex = v;
                    max_degree = degree;
                }
            }

            remaining_nodes.erase(std::find(remaining_nodes.begin(), remaining_nodes.end(), max_degree_vertex));

            std::vector<uint32_t> clique{max_degree_vertex};

            for (const uint32_t &u : remaining_nodes)
            {

                bool is_adjacent_to_all_vertices = true;
                for (const uint32_t &v : clique)
                {
                    if (!adjacency.areAdjacent(u, v))
                    {
                        is_adjacent_to_all_vertices = false;
                        break;
                    }
                }

                if (is_adjacent_to_all_vertices)
                {
                    clique.push_back(u);
                }
            }

            cliques.push_back(clique);
        }

        return getLargestListSize(cliques);
    }

    uint32_t findMaximalCliqueBronKerboschSimple()
    {
        std::vector<uint32_t> R;
        std::vector<uint32_t> P = nodes;
        std::vector<uint32_t> X;
        std::vector<std::vector<uint32_t>> maximalCliques;

        findMaximalCliqueBronKerboschSimple(R, P, X, maximalCliques);

        return getLargestListSize(maximalCliques);
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
                    if (adjacency.areAdjacent(P[i], P[j]))
                    {
                        newP.push_back(P[j]);
                    }
                }
                for (uint32_t j = 0; j < X.size(); j++)
                {
                    if (adjacency.areAdjacent(P[i], X[j]))
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

        return getLargestListSize(maximalCliques);
    }

    void findMaximalCliqueBronKerboschPivot(std::vector<uint32_t> &R, std::vector<uint32_t> &P, std::vector<uint32_t> &X, std::vector<std::vector<uint32_t>> &maximalCliques)
    {
        if (P.empty() && X.empty())
        {
            maximalCliques.push_back(R);
        }
        else
        {
            if (!P.empty())
            {
                uint32_t pivot = choosePivot(P, X);

                std::vector<uint32_t> P_without_neighbors_of_pivot;
                for (const uint32_t &v : P)
                {
                    if (!adjacency.areAdjacent(pivot, v))
                    {
                        P_without_neighbors_of_pivot.push_back(v);
                    }
                }

                for (const uint32_t &v : P_without_neighbors_of_pivot)
                {
                    std::vector<uint32_t> newR = R;
                    newR.push_back(v);
                    std::vector<uint32_t> newP;
                    std::vector<uint32_t> newX;
                    for (const uint32_t &u : P)
                    {
                        if (adjacency.areAdjacent(v, u))
                        {
                            newP.push_back(u);
                        }
                    }
                    for (const uint32_t &u : X)
                    {
                        if (adjacency.areAdjacent(v, u))
                        {
                            newX.push_back(u);
                        }
                    }
                    findMaximalCliqueBronKerboschPivot(newR, newP, newX, maximalCliques);
                    P.erase(std::remove(P.begin(), P.end(), v), P.end());
                    X.push_back(v);
                }
            }
        }
    }

    uint32_t choosePivot(const std::vector<uint32_t> &P, const std::vector<uint32_t> &X)
    {
        uint32_t pivot = P[0];
        uint32_t maxNeighbors = 0;
        for (const uint32_t &u : P)
        {
            uint32_t neighbors = 0;
            for (const uint32_t &v : P)
            {
                if (adjacency.areAdjacent(u, v))
                {
                    neighbors++;
                }
            }
            for (const uint32_t &v : X)
            {
                if (adjacency.areAdjacent(u, v))
                {
                    neighbors++;
                }
            }
            if (neighbors > maxNeighbors)
            {
                pivot = u;
                maxNeighbors = neighbors;
            }
        }
        return pivot;
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
        max = g.findMaximalCliqueGreedy();
        break;
    case 4:
        max = g.findMaximalCliqueBronKerboschSimple();
        break;
    case 5:
        max = g.findMaximalCliqueBronKerboschPivot();
        break;
    }

    timepoint_t end = clock_t::now();

    double elapsedSecondsAlgorithm = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0;
    std::cout << "Algorithm time: " << elapsedSecondsAlgorithm << "s" << std::endl;
    std::cout << "Maximum Clique Size: " << max << std::endl;
}

int main()
{

    std::cout << "Enter n and d to find the maximal clique of a graph with n nodes and d distance." << std::endl;
    std::cout << "Enter m for the method to use." << std::endl;
    std::cout << "Enter a to create adjacency matrix. y/n" << std::endl;
    std::cout << "Enter 1 for brute force." << std::endl;
    std::cout << "Enter 2 for Heuristic Breadth-First Search." << std::endl;
    std::cout << "Enter 3 for Heuristic Greedy." << std::endl;
    std::cout << "Enter 4 for Simple Bron-Kerbosch." << std::endl;
    std::cout << "Enter 5 for Pivot Bron-Kerbosch." << std::endl;
    std::cout << "Enter unvalid numbers to exit." << std::endl;

    while (true)
    {
        std::string temp;
        int n, d, m;

        std::cout << "Enter n: ";
        std::getline(std::cin, temp);
        if (!std::all_of(temp.begin(), temp.end(), [](char c)
                         { return std::isdigit(c); }) ||
            !(std::istringstream(temp) >> n) || n < 1)
        {
            break;
        }

        std::cout << "Enter d: ";
        std::getline(std::cin, temp);
        if (!std::all_of(temp.begin(), temp.end(), [](char c)
                         { return std::isdigit(c); }) ||
            !(std::istringstream(temp) >> d) || d < 1)
        {
            break;
        }

        std::cout << "Enter m: ";
        std::getline(std::cin, temp);
        if (!std::all_of(temp.begin(), temp.end(), [](char c)
                         { return std::isdigit(c); }) ||
            !(std::istringstream(temp) >> m) || m < 1 || m > 5)
        {
            break;
        }

        findMaximalClique(n, d, m);
    }

    return 0;
}

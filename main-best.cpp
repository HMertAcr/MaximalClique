#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <stack>
#include <tuple>
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
            adjacencies.push_back(std::unordered_set<uint32_t>());
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

    void findAdjacencies()
    {
        for (uint32_t i = 0; i < adjacencies.size(); i++)
        {
            for (uint32_t j = i + 1; j < adjacencies.size(); j++)
            {
                if (hammingDistance(i, j) >= distance)
                {
                    adjacencies[i].insert(j);
                    adjacencies[j].insert(i);
                }
            }
        }
        isAdjacencyCreated = true;
    }

    // Returns the set of d-adjacent nodes for the given node 'x'.
    const std::unordered_set<uint32_t> &GetAdjacencies(uint32_t x)
    {
        if (!isAdjacencyCreated)
        {
            if (adjacencies[x].size() == 0)
            {
                for (uint32_t i = 0; i < adjacencies.size(); i++)
                {
                    if (hammingDistance(x, i) >= distance)
                    {
                        adjacencies[x].insert(i);
                    }
                }
            }
        }
        return adjacencies.at(x);
    }

    // Returns true if nodes 'a' and 'b' are d-adjacent, false otherwise.
    bool areAdjacent(uint32_t a, uint32_t b) const
    {
        if (isAdjacencyCreated)
        {
            return adjacencies.at(a).count(b) > 0;
        }
        else
        {
            return hammingDistance(a, b) >= distance;
        }
        //could just calculate hamming distance even might be faster
    }

private:
    std::vector<std::unordered_set<uint32_t>> adjacencies;
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
        for (uint32_t i = 0; i < (1 << bitlength); i++)
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

    void createAdjacency()
    {
        adjacency.findAdjacencies();
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
        const std::unordered_set<uint32_t> &adjacencies = adjacency.GetAdjacencies(clique[0]);
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
                if (visited[i] == 0 && adjacency.areAdjacent(current, nodes[i]))
                {
                    bool is_connected_to_all = true;
                    for (const uint32_t &node : current_clique)
                    {
                        if (!adjacency.areAdjacent(current, nodes[i]))
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

    uint32_t findMaximalCliqueHeuristicDFS()
    {
        std::vector<std::vector<uint32_t>> cliques;
        std::vector<uint32_t> visited(nodes.size(), 0);
        std::vector<uint32_t> current_clique;

        for (uint32_t start = 0; start < nodes.size(); start++)
        {
            if (visited[start] == 0)
            {
                std::vector<uint32_t> stack = {start};
                visited[start] = 1;
                current_clique.push_back(start);
                while (!stack.empty())
                {
                    uint32_t current = stack.back();
                    stack.pop_back();
                    const std::unordered_set<uint32_t> &adjacencies = adjacency.GetAdjacencies(current);
                    for (const uint32_t &i : adjacencies)
                    {
                        if (visited[i] == 0)
                        {
                            bool is_connected_to_all = true;
                            for (const uint32_t &node : current_clique)
                            {
                                if (!adjacency.areAdjacent(nodes[i], node))
                                {
                                    is_connected_to_all = false;
                                    break;
                                }
                            }
                            if (is_connected_to_all)
                            {
                                visited[i] = 1;
                                current_clique.push_back(nodes[i]);
                                stack.push_back(i);
                            }
                        }
                    }
                }
                cliques.push_back(current_clique);
                current_clique.clear();
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

    uint32_t choosePivot(const std::vector<uint32_t> &P, const std::vector<uint32_t> &X)
    {
        uint32_t pivot = P[0];
        uint32_t minNeighbors = UINT32_MAX;
        for (uint32_t u : P)
        {
            uint32_t neighbors = 0;
            for (uint32_t v : P)
            {
                if (adjacency.areAdjacent(u, v))
                {
                    neighbors++;
                }
            }
            for (uint32_t v : X)
            {
                if (adjacency.areAdjacent(u, v))
                {
                    neighbors++;
                }
            }
            if (neighbors < minNeighbors)
            {
                pivot = u;
                minNeighbors = neighbors;
            }
        }
        return pivot;
    }

    void findMaximalCliqueBronKerboschPivot(std::vector<uint32_t> &R, std::vector<uint32_t> &P, std::vector<uint32_t> &X, std::vector<std::vector<uint32_t>> &maximalCliques)
    {
        if (P.empty() && X.empty())
        {
            maximalCliques.push_back(R);
        }
        else
        {
            uint32_t pivot = choosePivot(P, X);
            std::vector<uint32_t> P1;
            std::vector<uint32_t> X1;
            for (uint32_t v : P)
            {
                if (adjacency.areAdjacent(pivot, v))
                {
                    P1.push_back(v);
                }
            }
            for (uint32_t v : X)
            {
                if (adjacency.areAdjacent(pivot, v))
                {
                    X1.push_back(v);
                }
            }
            findMaximalCliqueBronKerboschPivot(R, P1, X1, maximalCliques);
            for (uint32_t u : P)
            {
                if (u != pivot)
                {
                    std::vector<uint32_t> newR = R;
                    newR.push_back(u);
                    std::vector<uint32_t> newP;
                    std::vector<uint32_t> newX;
                    for (uint32_t v : P)
                    {
                        if (adjacency.areAdjacent(u, v))
                        {
                            newP.push_back(v);
                        }
                    }
                    for (uint32_t v : X)
                    {
                        if (adjacency.areAdjacent(u, v))
                        {
                            newX.push_back(v);
                        }
                    }
                    findMaximalCliqueBronKerboschPivot(newR, newP, newX, maximalCliques);
                    P.erase(std::remove(P.begin(), P.end(), u), P.end());
                    X.push_back(u);
                }
            }
        }
    }

    uint32_t findMaximalCliqueTTT()
    {
        std::vector<std::unordered_set<uint32_t>> cliques;
        std::vector<uint32_t> remaining_nodes = nodes;

        while (!remaining_nodes.empty())
        {
            // Choose a vertex with the maximum degree.
            uint32_t max_degree_vertex = remaining_nodes[0];

            uint32_t max_degree = adjacency.GetAdjacencies(max_degree_vertex).size();

            for (uint32_t v : remaining_nodes)
            {
                uint32_t degree = adjacency.GetAdjacencies(v).size();
                if (degree > max_degree)
                {
                    max_degree_vertex = v;
                    max_degree = degree;
                }
            }

            remaining_nodes.erase(std::find(remaining_nodes.begin(), remaining_nodes.end(), max_degree_vertex));

            // Initialize the clique to be {v}.
            std::unordered_set<uint32_t> clique{max_degree_vertex};

            // Add all vertices that are adjacent to all vertices in the clique.
            for (uint32_t u : remaining_nodes)
            {

                bool is_adjacent_to_all_vertices = true;
                for (uint32_t v : clique)
                {
                    if (!adjacency.areAdjacent(u, v))
                    {
                        is_adjacent_to_all_vertices = false;
                        break;
                    }
                }

                if (is_adjacent_to_all_vertices)
                {
                    clique.insert(u);
                }
            }

            // Add the clique to the set of cliques.
            cliques.push_back(clique);
        }

        uint32_t max = 0;
        for (uint32_t i = 0; i < cliques.size(); i++)
        {
            if (cliques[i].size() > max)
            {
                max = cliques[i].size();
            }
        }

        return max;
    }
};

void findMaximalClique(const uint32_t n, const uint32_t d, const uint32_t m, std::string a)
{

    using clock_t = std::chrono::high_resolution_clock;
    using timepoint_t = std::chrono::time_point<clock_t>;

    timepoint_t start = clock_t::now();

    graph g(n, d);

    if (a == "y" || a == "Y")
    {
        g.createAdjacency();
    }
    timepoint_t midpoint = clock_t::now();

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
        max = g.findMaximalCliqueHeuristicDFS();
        break;
    case 4:
        max = g.findMaximalCliqueBronKerboschSimple();
        break;
    case 5:
        max = g.findMaximalCliqueBronKerboschPivot();
        break;
    case 6:
        max = g.findMaximalCliqueTTT();
        break;
    }

    timepoint_t end = clock_t::now();

    std::cout << "Maximum Clique Size: " << max << std::endl;

    if (a == "y" || a == "Y")
    {
        double elapsedSecondsGraphAndAdjacency = std::chrono::duration_cast<std::chrono::microseconds>(midpoint - start).count() / 1000000.0;
        std::cout << "Time spent creating graph and adjacencies: " << elapsedSecondsGraphAndAdjacency << "s" << std::endl;
    }

    double elapsedSecondsAlgorithm = std::chrono::duration_cast<std::chrono::microseconds>(end - midpoint).count() / 1000000.0;
    std::cout << "Algorithm time: " << elapsedSecondsAlgorithm << "s" << std::endl;
}

int main()
{

    std::cout << "Enter n and d to find the maximal clique of a graph with n nodes and d distance." << std::endl;
    std::cout << "Enter m for the method to use." << std::endl;
    std::cout << "Enter a to create adjacency matrix. y/n" << std::endl;
    std::cout << "Enter 1 for brute force." << std::endl;
    std::cout << "Enter 2 for Heuristic Breadth-First Search." << std::endl;
    std::cout << "Enter 3 for Heuristic Depth-First Search." << std::endl;
    std::cout << "Enter 4 for Simple Bron-Kerbosch." << std::endl;
    std::cout << "Enter 5 for Pivot Bron-Kerbosch." << std::endl;
    std::cout << "Enter 6 for Tomita, Tanaka, and Takahashi (TTT)." << std::endl;
    std::cout << "Enter unvalid numbers to exit." << std::endl;

    while (true)
    {

        std::string temp;

        std::cout << "Enter n: ";
        std::getline(std::cin, temp);

        bool isNumber = true;
        for (char c : temp)
        {
            if (!isdigit(c))
            {
                isNumber = false;
                break;
            }
        }

        if (!isNumber || stoi(temp) < 1)
        {
            break;
        }
        uint32_t n = stoi(temp);

        std::cout << "Enter d: ";
        std::getline(std::cin, temp);

        for (char c : temp)
        {
            if (!isdigit(c))
            {
                isNumber = false;
                break;
            }
        }

        if (!isNumber || stoi(temp) < 1)
        {
            break;
        }
        uint32_t d = stoi(temp);

        std::cout << "Enter m: ";
        std::getline(std::cin, temp);

        for (char c : temp)
        {
            if (!isdigit(c))
            {
                isNumber = false;
                break;
            }
        }

        if (!isNumber || stoi(temp) < 1 || stoi(temp) > 6)
        {
            break;
        }
        uint32_t m = stoi(temp);

        std::cout << "Enter a: ";
        std::getline(std::cin, temp);

        if (!(temp == "y" || temp == "Y" || temp == "n" || temp == "N"))
        {
            break;
        }

        std::string a = temp;

        findMaximalClique(n, d, m, a);
    }

    return 0;
}
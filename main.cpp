#include <iostream>
#include <string>
#include <vector>

std::string binaryRep(const unsigned int bit, const unsigned int bitlength)
{
    std::string s;
    s.reserve(bitlength);
    s = "";

    for (unsigned int i = 0; i < bitlength; i++)
    {
        s = ((bit & (1 << i)) ? "1" : "0") + s;
    }
    return s;
}

static unsigned int hammingDistance(const unsigned int bit_a, const unsigned int bit_b)
{
    int diff = bit_a ^ bit_b;
    int count = 0;

    while (diff != 0)
    {
        count++;
        diff &= diff - 1;
    }

    return count;
}

class graph
{
public:
    std::vector<unsigned int> nodes;
    unsigned int bitlength;

    // Default constructor
    graph()
    {
    }

    // Constructor
    graph(unsigned int bitlength)
    {
        this->bitlength = bitlength;

        for (unsigned int i = 0; i < (1 << bitlength); i++)
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

    static bool isSameNodes(std::vector<unsigned int> nodes_a, std::vector<unsigned int> nodes_b)
    {
        if (nodes_a.size() != nodes_b.size())
        {
            return false;
        }

        for (unsigned int i = 0; i < nodes_a.size(); i++)
        {
            bool found = false;
            for (unsigned int j = 0; j < nodes_b.size(); j++)
            {
                if (nodes_a[i] == nodes_b[j])
                {
                    nodes_b.erase(nodes_b.begin() + j);
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

    static void removeDuplicateNodes(std::vector<std::vector<unsigned int>> &nodes)
    {

        std::vector<std::vector<unsigned int>> newNodes;

        for (unsigned int i = 0; i < nodes.size(); i++)
        {
            bool found = false;
            for (unsigned int j = 0; j < newNodes.size(); j++)
            {
                if (isSameNodes(nodes[i], newNodes[j]))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                newNodes.push_back(nodes[i]);
            }
        }

        nodes = newNodes;
    }

    std::vector<std::vector<unsigned int>> getSubsets(std::vector<unsigned int> &nums)
    {

        std::vector<std::vector<unsigned int>> subset;
        std::vector<unsigned int> empty;
        subset.push_back(empty);

        for (int i = 0; i < nums.size(); i++)
        {
            std::vector<std::vector<unsigned int>> subsetTemp = subset;

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

    bool checkClique(std::vector<unsigned int> clique, unsigned int d)
    {
        for (unsigned int i = 0; i < clique.size(); i++)
        {
            for (unsigned int j = i + 1; j < clique.size(); j++)
            {
                if (hammingDistance(clique[i], clique[j]) < d)
                {
                    return false;
                }
            }
        }
        return true;
    }

    bool checkMaxmialClique(std::vector<unsigned int> clique, unsigned int d)
    {
        for (unsigned int i = 0; i < nodes.size(); i++)
        {
            bool found = false;
            for (unsigned int j = 0; j < clique.size(); j++)
            {
                if (nodes[i] == clique[j])
                {
                    found = true;
                    break;
                }
            }
            if (!found)
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

    void findMaximalCliqueWorstBF(const unsigned int d)
    {
        // Generate all possible subsets of nodes
        std::vector<std::vector<unsigned int>> subsets = getSubsets(nodes);

        // Remove all subsets that are not cliques
        for (unsigned int i = 0; i < subsets.size(); i++)
        {
            if (!checkClique(subsets[i], d))
            {
                subsets.erase(subsets.begin() + i);
                i--;
            }
        }

        // Remove all subsets that are not maximal cliques
        for (unsigned int i = 0; i < subsets.size(); i++)
        {
            if (!checkMaxmialClique(subsets[i], d))
            {
                subsets.erase(subsets.begin() + i);
                i--;
            }
        }

        unsigned int max = 0;
        for (unsigned int i = 0; i < subsets.size(); i++)
        {
            if (subsets[i].size() > max)
            {
                max = subsets[i].size();
            }
        }

        // unsigned int count = 0;
        // for (unsigned int i = 0; i < subsets.size(); i++)
        // {
        //     if (subsets[i].size() == max)
        //     {
        //         count++;
        //     }
        // }

        // std::cout << "Maximal Clique Count: " << subsets.size() << std::endl;
        // std::cout << "Maximum Clique Count: " << count << std::endl;

        std::cout << "Maximum Clique Size: " << max << std::endl;

    }

    void findMaximalCliqueBetterBF(const unsigned int d)
    {
        std::vector<std::vector<unsigned int>> clique;

        // Find all Cliques
        for (unsigned int i = 0; i < nodes.size(); i++)
        {
            clique.push_back(std::vector<unsigned int>());
            clique[i].push_back(nodes[i]);
            for (unsigned int j = 0; j < nodes.size(); j++)
            {
                bool canAdd = true;
                for (unsigned int k = 0; k < clique[i].size(); k++)
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

        graph::removeDuplicateNodes(clique);

        unsigned int max = 0;
        for (unsigned int i = 0; i < clique.size(); i++)
        {
            if (clique[i].size() > max)
            {
                max = clique[i].size();
            }
        }

        // unsigned int count = 0;
        // for (unsigned int i = 0; i < clique.size(); i++)
        // {
        //     if (clique[i].size() == max)
        //     {
        //         count++;
        //     }
        // }

        // std::cout << "Maximal Clique Count: " << clique.size() << std::endl;
        // std::cout << "Maximum Clique Count: " << count << std::endl;

        std::cout << "Maximum Clique Size: " << max << std::endl;
        
    }
};

void findMaximalClique(const unsigned int n, const unsigned int d, const unsigned int m)
{
    graph g(n);

    if (m == 0)
    {
        g.findMaximalCliqueWorstBF(d);
    }
    else if (m == 1)
    {
        g.findMaximalCliqueBetterBF(d);
    }
}

int main()
{

    std::cout << "Enter n and d to find the maximal clique of a graph with n nodes and d distance." << std::endl;
    std::cout << "Enter m for the method to use." << std::endl;
    std::cout << "Enter 0 for worst brute force." << std::endl;
    std::cout << "Enter 1 for better brute force." << std::endl;
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
        unsigned int n = stoi(temp);

        std::cout << "Enter d: ";
        std::getline(std::cin, temp);

        if (stoi(temp) < 1)
        {
            break;
        }
        unsigned int d = stoi(temp);

        std::cout << "Enter m: ";
        std::getline(std::cin, temp);

        if (stoi(temp) < 0 || stoi(temp) > 1)
        {
            break;
        }
        unsigned int m = stoi(temp);

        findMaximalClique(n, d, m);
    }

    return 0;
}
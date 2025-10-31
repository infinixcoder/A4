#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <ctime>

using namespace std;

// --- Data Structures ---

class BayesNode
{
private:
    string nodeName;
    vector<int> childIndices;
    vector<int> parentIndices;
    int valueCount;
    vector<string> valueNames;
    map<string, int> valueIndexMap;
    vector<float> cptTable;

public:
    BayesNode(string name, int n, vector<string> vals)
    {
        nodeName = name;
        valueCount = n;
        valueNames = vals;
        for (int i = 0; i < valueCount; ++i)
        {
            valueIndexMap[valueNames[i]] = i;
        }
    }

    string getName() const
    {
        return nodeName;
    }

    vector<int> getChildren() const
    {
        return childIndices;
    }

    vector<int> getParentIndices() const
    {
        return parentIndices;
    }

    vector<float> &getCPT()
    {
        return cptTable;
    }

    const vector<float> &getCPT() const
    {
        return cptTable;
    }

    int getValueCount() const
    {
        return valueCount;
    }

    vector<string> getValues() const
    {
        return valueNames;
    }

    int findValueIndex(const string &val) const
    {
        try
        {
            return valueIndexMap.at(val);
        }
        catch (const out_of_range &e)
        {
            if (val.length() >= 2 && val.front() == '"' && val.back() == '"')
            {
                string stripped_val = val.substr(1, val.length() - 2);
                try
                {
                    return valueIndexMap.at(stripped_val);
                }
                catch (const out_of_range &e2)
                {
                }
            }
            cerr << "Error: Value '" << val << "' not found in node '" << nodeName << "'" << endl;
            return -1;
        }
    }

    void setCPT(vector<float> newCPT)
    {
        cptTable.clear();
        cptTable = newCPT;
    }

    void setParentIndices(vector<int> newParentIndices)
    {
        parentIndices.clear();
        parentIndices = newParentIndices;
    }

    int addChild(int newChildIndex)
    {
        for (int i = 0; i < (int)childIndices.size(); i++)
        {
            if (childIndices[i] == newChildIndex)
                return 0;
        }
        childIndices.push_back(newChildIndex);
        return 1;
    }
};

class BayesNetwork
{
public:
    vector<BayesNode> allNodes;
    map<string, int> nameToIndexMap;

public:
    int addNode(BayesNode node)
    {
        int index = allNodes.size();
        allNodes.push_back(node);
        nameToIndexMap[node.getName()] = index;
        return 0;
    }

    list<BayesNode>::iterator getNode(int i)
    {
        static list<BayesNode> dummy_list;
        return dummy_list.begin();
    }

    BayesNode &getNodeByIndex(int i)
    {
        return allNodes[i];
    }

    const BayesNode &getNodeByIndex(int i) const
    {
        return allNodes[i];
    }

    BayesNode &getNodeByName(string val_name)
    {
        return allNodes[nameToIndexMap.at(val_name)];
    }

    int getNetSize() const
    {
        return allNodes.size();
    }

    int getIndexForName(string val_name) const
    {
        if (nameToIndexMap.count(val_name))
        {
            return nameToIndexMap.at(val_name);
        }
        return -1;
    }
};

// --- Core Logic Helpers ---

string trimWhitespace(const string &str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

int findFlatCptIndex(const BayesNetwork &bn, int nodeIndex, const vector<string> &dataRow)
{
    const BayesNode &node = bn.getNodeByIndex(nodeIndex);
    int nValues = node.getValueCount();

    int valueIdx = node.findValueIndex(dataRow[nodeIndex]);
    if (valueIdx == -1)
        return -1;

    vector<int> parentIdxs = node.getParentIndices();
    if (parentIdxs.empty())
    {
        return valueIdx;
    }

    int parentConfigIndex = 0;
    int multiplier = 1;

    for (int i = (int)parentIdxs.size() - 1; i >= 0; --i)
    {
        int pIndex = parentIdxs[i];
        const BayesNode &parentNode = bn.getNodeByIndex(pIndex);

        string pValue = dataRow[pIndex];
        int pValueIndex = parentNode.findValueIndex(pValue);
        if (pValueIndex == -1)
            return -1;

        parentConfigIndex += pValueIndex * multiplier;
        multiplier *= parentNode.getValueCount();
    }

    int finalIndex = (parentConfigIndex * nValues) + valueIdx;
    return finalIndex;
}

float computeRecordJointProb(const BayesNetwork &bn, const vector<string> &dataRow)
{
    double logProbSum = 0.0;
    for (int i = 0; i < bn.getNetSize(); ++i)
    {
        int cptIdx = findFlatCptIndex(bn, i, dataRow);
        if (cptIdx == -1)
        {
            cerr << "Error calculating joint prob: invalid dataRow." << endl;
            return 0.0;
        }

        float prob = bn.getNodeByIndex(i).getCPT()[cptIdx];
        if (prob == 0.0f)
        {
            return 0.0;
        }
        logProbSum += log(prob);
    }

    return exp(logProbSum);
}

// --- EM Algorithm Steps ---

// FIX 1: Changed to return map<int, vector<double>>
map<int, vector<double>> performEStep(const BayesNetwork &bn, const vector<vector<string>> &dataEntries)
{
    // FIX 1: Use vector<double> for high-precision accumulation
    map<int, vector<double>> sufficientStatistics;
    for (int i = 0; i < bn.getNetSize(); ++i)
    {
        int cptSize = bn.getNodeByIndex(i).getCPT().size();
        // FIX 1: Use vector<double> and initialize with double literal 0.6
        sufficientStatistics[i] = vector<double>(cptSize, 0.6);
    }

    for (const auto &row : dataEntries)
    {
        int unknownVarIndex = 0;
        while (unknownVarIndex < (int)row.size() && row[unknownVarIndex] != "?")
        {
            unknownVarIndex++;
        }

        if (unknownVarIndex == (int)row.size())
        {
            for (int i = 0; i < bn.getNetSize(); ++i)
            {
                int cptIndex = findFlatCptIndex(bn, i, row);
                if (cptIndex != -1)
                {
                    // FIX 1: Add 1.0 (double) to vector<double>
                    sufficientStatistics[i][cptIndex] += 1.0;
                }
            }
        }
        else
        {
            const BayesNode &missingNode = bn.getNodeByIndex(unknownVarIndex);
            vector<string> possibleValues = missingNode.getValues();
            vector<float> posterior(possibleValues.size());
            float totalProb = 0.0f;

            for (int i = 0; i < (int)possibleValues.size(); ++i)
            {
                vector<string> tempRow = row;
                tempRow[unknownVarIndex] = possibleValues[i];
                float jointProb = computeRecordJointProb(bn, tempRow);
                posterior[i] = jointProb;
                totalProb += jointProb;
            }

            if (totalProb > 0)
            {
                for (int i = 0; i < (int)posterior.size(); ++i)
                {
                    posterior[i] /= totalProb;
                }
            }
            else
            {
                for (int i = 0; i < (int)posterior.size(); ++i)
                {
                    posterior[i] = 1.0f / posterior.size();
                }
            }

            for (int i = 0; i < (int)possibleValues.size(); ++i)
            {
                float weight = posterior[i];
                if (weight == 0.0f)
                {
                    continue;
                }

                vector<string> tempRow = row;
                tempRow[unknownVarIndex] = possibleValues[i];

                for (int j = 0; j < bn.getNetSize(); ++j)
                {
                    int cptIndex = findFlatCptIndex(bn, j, tempRow);
                    if (cptIndex != -1)
                    {
                        // FIX 1: Add float 'weight' to double accumulator
                        sufficientStatistics[j][cptIndex] += weight;
                    }
                }
            }
        }
    }
    return sufficientStatistics;
}

// FIX 1: Changed parameter to map<int, vector<double>>
float performMStep(BayesNetwork &bn, const map<int, vector<double>> &sufficientStatistics)
{
    float maxDelta = 0.0f;

    for (int i = 0; i < bn.getNetSize(); ++i)
    {
        BayesNode &node = bn.getNodeByIndex(i);
        vector<float> &cptTable = node.getCPT();
        // FIX 1: Use const vector<double>&
        const vector<double> &counts = sufficientStatistics.at(i);
        int valueCount = node.getValueCount();
        int numParentConfigs = cptTable.size() / valueCount;

        for (int configIdx = 0; configIdx < numParentConfigs; ++configIdx)
        {
            int baseIdx = configIdx * valueCount;
            // FIX 1: Use double for sum
            double sumOfCounts = 0.0;
            for (int k = 0; k < valueCount; ++k)
            {
                sumOfCounts += counts[baseIdx + k];
            }

            for (int k = 0; k < valueCount; ++k)
            {
                // FIX 1: Use double for newProb
                double newProb;
                if (sumOfCounts == 0.0)
                {
                    // FIX 1: Use double literal 1.0
                    newProb = 1.0 / valueCount;
                }
                else
                {
                    newProb = counts[baseIdx + k] / sumOfCounts;
                }

                // Change is calculated as (double - float) -> double
                float change = abs(newProb - cptTable[baseIdx + k]);
                if (change > maxDelta)
                {
                    maxDelta = change;
                }
                // Store final result as float
                cptTable[baseIdx + k] = newProb;
            }
        }
    }
    return maxDelta;
}

void executeEMAlgorithm(BayesNetwork &bn, const vector<vector<string>> &dataEntries, int iterLimit, float epsilon)
{
    cout << "Starting EM algorithm... (Refactored)" << endl;
    cout << "Max Iterations: " << iterLimit << ", Convergence: " << epsilon << endl;

    for (int iter = 0; iter < iterLimit; ++iter)
    {
        // FIX 1: Type of sufficientStatistics is now map<int, vector<double>>
        map<int, vector<double>> sufficientStatistics = performEStep(bn, dataEntries);
        float maxDelta = performMStep(bn, sufficientStatistics);

        cout << "Iteration " << iter + 1 << ": Max param change = " << maxDelta << endl;
        if (maxDelta < epsilon)
        {
            cout << "Convergence reached!" << endl;
            break;
        }
    }
    cout << "EM algorithm finished." << endl;
}

// --- File I/O Functions ---

BayesNetwork loadNetworkFromFile(const char *filePath)
{
    BayesNetwork bn;
    string line;
    ifstream netFile(filePath);

    if (!netFile.is_open())
    {
        cout << "Error: Could not open file " << filePath << endl;
        return bn;
    }

    while (getline(netFile, line))
    {
        line = trimWhitespace(line);
        if (line.empty() || line[0] == '#')
            continue;

        stringstream ss(line);
        string token;
        ss >> token;

        if (token == "variable")
        {
            string varName;
            ss >> varName;
            getline(netFile, line);
            stringstream ss2(line);
            string k1, k2, k3, k4, k5, k6;
            int numValues;
            ss2 >> k1 >> k2 >> k3 >> numValues >> k4 >> k5 >> k6;

            vector<string> values;
            string value;
            while (ss2 >> value)
            {
                if (value == "};")
                    break;
                if (value.back() == ',')
                {
                    value = value.substr(0, value.length() - 1);
                }
                values.push_back(value);
            }
            BayesNode newNode(varName, numValues, values);
            bn.addNode(newNode);
        }
        else if (token == "probability")
        {
            string paren, nodeName;
            ss >> paren >> nodeName;

            string fullProbLine = line;
            while (fullProbLine.find('{') == string::npos)
            {
                string nextLine;
                if (!getline(netFile, nextLine))
                    break;
                fullProbLine += " " + trimWhitespace(nextLine);
            }

            size_t startParen = fullProbLine.find('(');
            size_t endParen = fullProbLine.find(')');
            size_t pipePos = fullProbLine.find('|');

            if (startParen == string::npos || endParen == string::npos)
                continue;

            string probContent = fullProbLine.substr(startParen + 1, endParen - startParen - 1);
            stringstream probSS(probContent);
            probSS >> nodeName;

            int index = bn.getIndexForName(nodeName);
            BayesNode &node = bn.getNodeByIndex(index);

            vector<int> parentIndices;
            if (pipePos != string::npos && pipePos < endParen)
            {
                string parentsStr = fullProbLine.substr(pipePos + 1, endParen - pipePos - 1);
                stringstream parentSS(parentsStr);
                string parentName;
                while (parentSS >> parentName)
                {
                    if (parentName.back() == ',')
                    {
                        parentName = parentName.substr(0, parentName.length() - 1);
                    }
                    int parentIndex = bn.getIndexForName(parentName);
                    parentIndices.push_back(parentIndex);
                    bn.getNodeByIndex(parentIndex).addChild(index);
                }
            }
            node.setParentIndices(parentIndices);

            vector<float> cptData;
            while (getline(netFile, line))
            {
                line = trimWhitespace(line);
                if (line == "};")
                    break;
                if (line.empty())
                    continue;

                size_t closeParen = line.find(')');
                string probPart;
                if (closeParen != string::npos)
                {
                    probPart = line.substr(closeParen + 1);
                }
                else if (line.find("table") != string::npos)
                {
                    probPart = line.substr(line.find("table") + 5);
                }
                else
                {
                    probPart = line;
                }

                stringstream ssProb(probPart);
                string numToken;
                while (ssProb >> numToken)
                {
                    while (!numToken.empty() && (numToken.back() == ',' || numToken.back() == ';'))
                    {
                        numToken = numToken.substr(0, numToken.length() - 1);
                    }
                    if (!numToken.empty() && (isdigit(numToken[0]) || numToken[0] == '.' || numToken[0] == '-'))
                    {
                        cptData.push_back(atof(numToken.c_str()));
                    }
                }
            }
            node.setCPT(cptData);
        }
    }
    netFile.close();
    return bn;
}

void saveNetworkToFile(const char *filePath, BayesNetwork &bn)
{
    ofstream outFile(filePath);
    if (!outFile.is_open())
    {
        cout << "Error: Could not open file " << filePath << " for writing" << endl;
        return;
    }

    outFile << "// Bayesian Network (Learned by EM)" << endl
            << endl;

    int N = bn.getNetSize();
    for (int i = 0; i < N; i++)
    {
        BayesNode &node = bn.getNodeByIndex(i);
        outFile << "variable " << node.getName() << " {" << endl;
        outFile << "  type discrete [ " << node.getValueCount() << " ] = { ";
        vector<string> vals = node.getValues();
        for (int j = 0; j < (int)vals.size(); j++)
        {
            outFile << vals[j] << (j < (int)vals.size() - 1 ? ", " : "");
        }
        outFile << " };" << endl;
        outFile << "}" << endl;
    }

    outFile << std::fixed << std::setprecision(4);
    for (int i = 0; i < N; i++)
    {
        BayesNode &node = bn.getNodeByIndex(i);
        vector<int> parentIndices = node.getParentIndices();
        vector<string> values = node.getValues();
        vector<float> cpt = node.getCPT();

        outFile << "probability ( " << node.getName();
        if (!parentIndices.empty())
        {
            outFile << " | ";
            for (int j = 0; j < (int)parentIndices.size(); j++)
            {
                outFile << bn.getNodeByIndex(parentIndices[j]).getName() << (j < (int)parentIndices.size() - 1 ? ", " : "");
            }
        }
        outFile << " ) {" << endl;

        vector<int> radices;
        radices.reserve(parentIndices.size());
        for (int pIndex : parentIndices)
        {
            radices.push_back(bn.getNodeByIndex(pIndex).getValueCount());
        }

        int parentCombinations = 1;
        for (int r : radices)
            parentCombinations *= r;

        int cptIndex = 0;
        if (parentIndices.empty())
        {
            outFile << "    table ";
            for (int k = 0; k < (int)values.size(); k++)
            {
                outFile << (cptIndex < (int)cpt.size() ? cpt[cptIndex++] : -1) << (k < (int)values.size() - 1 ? ", " : "");
            }
            outFile << ";" << endl;
        }
        else
        {
            for (int comb = 0; comb < parentCombinations; comb++)
            {
                vector<int> idx(parentIndices.size(), 0);
                int tmp = comb;
                for (int p = (int)parentIndices.size() - 1; p >= 0; p--)
                {
                    idx[p] = tmp % radices[p];
                    tmp /= radices[p];
                }

                outFile << "    ( ";
                for (int p = 0; p < (int)parentIndices.size(); p++)
                {
                    BayesNode &pnode = bn.getNodeByIndex(parentIndices[p]);
                    outFile << pnode.getValues()[idx[p]] << (p < (int)parentIndices.size() - 1 ? ", " : "");
                }
                outFile << " ) ";

                for (int k = 0; k < (int)values.size(); k++)
                {
                    outFile << (cptIndex < (int)cpt.size() ? cpt[cptIndex++] : -1) << (k < (int)values.size() - 1 ? ", " : "");
                }
                outFile << ";" << endl;
            }
        }
        outFile << "};" << endl
                << endl;
    }
    outFile.close();
}

vector<vector<string>> loadDataRecords(const char *filePath)
{
    vector<vector<string>> dataEntries;
    ifstream file(filePath);
    string line;

    if (!file.is_open())
    {
        cout << "Error: Could not open data file " << filePath << endl;
        return dataEntries;
    }

    while (getline(file, line))
    {
        line = trimWhitespace(line);
        if (line.empty())
            continue;

        vector<string> row;
        stringstream ss(line);
        string value;

        while (getline(ss, value, ','))
        {
            if (!value.empty() && value.front() == '"')
            {
                value = value.substr(1);
            }
            if (!value.empty() && value.back() == '"')
            {
                value = value.substr(0, value.length() - 1);
            }
            row.push_back(value);
        }

        if (row.size() == 56)
        {
            dataEntries.push_back(row);
        }
        else if (!row.empty())
        {
            cout << "Warning: Skipping malformed row with " << row.size() << " fields, expected 56." << endl;
        }
    }
    file.close();
    return dataEntries;
}

// FIX 2: Replaced the entire function with the corrected logic
void initRandomProbs(BayesNetwork &bn)
{
    srand(time(0)); // Seed random number generator
    for (int i = 0; i < bn.getNetSize(); ++i)
    {
        BayesNode &node = bn.getNodeByIndex(i);
        vector<float> &cptTable = node.getCPT();
        int nValues = node.getValueCount();

        // 1. Calculate the *expected* CPT size
        int parentCombos = 1;
        for (int pIdx : node.getParentIndices())
        {
            parentCombos *= bn.getNodeByIndex(pIdx).getValueCount();
        }
        int expectedSize = parentCombos * nValues;

        // 2. Check if the CPT is just a placeholder (e.g., {-1})
        bool isPlaceholder = false;
        if (!cptTable.empty() && cptTable[0] == -1.0f)
        {
            isPlaceholder = true;
        }

        // 3. If it's empty OR a placeholder with the wrong size,
        //    reset the vector to the correct size and fill with -1.
        if (cptTable.empty() || (isPlaceholder && cptTable.size() != expectedSize))
        {
            // Use assign() to completely reset the vector's size and content
            cptTable.assign(expectedSize, -1.0f);
        }

        // 4. Now, this loop is safe because cptTable is guaranteed
        //    to have the correct size (expectedSize).
        for (int j = 0; j < (int)cptTable.size(); j += nValues)
        {
            // We only fill blocks that are marked as unknown.
            if (cptTable[j] == -1.0f)
            {
                float sum = 0.0f;
                vector<float> randomVals(nValues);
                for (int k = 0; k < nValues; ++k)
                {
                    randomVals[k] = (float)rand() / RAND_MAX + 0.001f;
                    sum += randomVals[k];
                }
                // Normalize and fill the CPT block
                for (int k = 0; k < nValues; ++k)
                {
                    cptTable[j + k] = randomVals[k] / sum;
                }
            }
        }
    }
}

// --- Main Driver ---

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cout << "Usage: ./run.sh <bif_file> <data_file>" << endl;
        return 1;
    }

    string bifFilePath = argv[1];
    string dataFilePath = argv[2];
    string outputFilePath = "solved_hailfinder.bif";

    BayesNetwork bayesNet = loadNetworkFromFile(bifFilePath.c_str());
    if (bayesNet.getNetSize() == 0)
    {
        cout << "Error: Failed to read network." << endl;
        return 1;
    }

    vector<vector<string>> dataRecords = loadDataRecords(dataFilePath.c_str());
    if (dataRecords.empty())
    {
        cout << "Error: Failed to read data or data file is empty." << endl;
        return 1;
    }

    initRandomProbs(bayesNet);

    executeEMAlgorithm(bayesNet, dataRecords, 50, 2e-7);

    saveNetworkToFile(outputFilePath.c_str(), bayesNet);

    return 0;
}
#include <iostream>
#include <string>
#include <vector>
#include <list> // Still included for parser compatibility, but core structure is vector
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <cmath>
#include <stdexcept> // For .at()
#include <ctime>     // For random seed

using namespace std;

// --- Data Structures (Refactored for Performance) ---

class Graph_Node
{
private:
    string Node_Name;
    vector<int> Children;
    vector<int> Parent_Indices; // Changed from vector<string>
    int nvalues;
    vector<string> values;
    map<string, int> value_to_index; // Added for O(1) value index lookup
    vector<float> CPT;

public:
    Graph_Node(string name, int n, vector<string> vals)
    {
        Node_Name = name;
        nvalues = n;
        values = vals;
        // Populate the value_to_index map
        for (int i = 0; i < nvalues; ++i)
        {
            value_to_index[values[i]] = i;
        }
    }

    string get_name() const
    {
        return Node_Name;
    }

    vector<int> get_children() const
    {
        return Children;
    }

    vector<int> get_Parent_Indices() const
    { // Changed
        return Parent_Indices;
    }

    vector<float> &get_CPT()
    { // Return by ref for direct modification
        return CPT;
    }

    // Const version for non-mutating access
    const vector<float> &get_CPT() const
    {
        return CPT;
    }

    int get_nvalues() const
    {
        return nvalues;
    }

    vector<string> get_values() const
    {
        return values;
    }

    int get_value_index(const string &val) const
    {
        try
        {
            return value_to_index.at(val);
        }
        catch (const out_of_range &e)
        {
            // A common issue is extra quotes, let's try to strip them
            if (val.length() >= 2 && val.front() == '"' && val.back() == '"')
            {
                string stripped_val = val.substr(1, val.length() - 2);
                try
                {
                    return value_to_index.at(stripped_val);
                }
                catch (const out_of_range &e2)
                {
                }
            }
            cerr << "Error: Value '" << val << "' not found in node '" << Node_Name << "'" << endl;
            return -1; // Or throw
        }
    }

    void set_CPT(vector<float> new_CPT)
    {
        CPT.clear();
        CPT = new_CPT;
    }

    void set_Parent_Indices(vector<int> Parent_Node_Indices)
    { // Changed
        Parent_Indices.clear();
        Parent_Indices = Parent_Node_Indices;
    }

    int add_child(int new_child_index)
    {
        for (int i = 0; i < (int)Children.size(); i++)
        {
            if (Children[i] == new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }
};

class network
{
public:
    vector<Graph_Node> nodes;       // Changed from list to vector for O(1) indexing
    map<string, int> name_to_index; // Added for O(log n) name lookup

public:
    int addNode(Graph_Node node)
    {
        int index = nodes.size();
        nodes.push_back(node);
        name_to_index[node.get_name()] = index;
        return 0;
    }

    // This function is just a stub to make the old parser compile.
    // It's not used in the EM logic.
    list<Graph_Node>::iterator getNode(int i)
    {
        static list<Graph_Node> dummy_list;
        return dummy_list.begin();
    }

    Graph_Node &get_node_by_index(int i)
    {
        return nodes[i];
    }

    const Graph_Node &get_node_by_index(int i) const
    {
        return nodes[i];
    }

    Graph_Node &get_node_by_name(string val_name)
    {
        return nodes[name_to_index.at(val_name)];
    }

    int netSize() const
    {
        return nodes.size();
    }

    int get_index(string val_name) const
    {
        if (name_to_index.count(val_name))
        {
            return name_to_index.at(val_name);
        }
        return -1;
    }
};

// --- Parser (Refactored to use new network class) ---

string trim(const string &str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

network read_network(const char *filename)
{
    network BayesNet;
    string line;
    ifstream myfile(filename);

    if (!myfile.is_open())
    {
        cout << "Error: Could not open file " << filename << endl;
        return BayesNet;
    }

    while (getline(myfile, line))
    {
        line = trim(line);
        if (line.empty() || line[0] == '#')
            continue;

        stringstream ss(line);
        string token;
        ss >> token;

        if (token == "variable")
        {
            string var_name;
            ss >> var_name;

            getline(myfile, line);
            stringstream ss2(line);

            string type_keyword, discrete_keyword, bracket, equals;
            int num_values;
            ss2 >> type_keyword >> discrete_keyword >> bracket >> num_values >> bracket >> equals >> bracket;

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
            Graph_Node new_node(var_name, num_values, values);
            BayesNet.addNode(new_node);
        }
        else if (token == "probability")
        {
            string paren, node_name;
            ss >> paren >> node_name;

            string full_prob_line = line;
            while (full_prob_line.find('{') == string::npos)
            {
                string next_line;
                if (!getline(myfile, next_line))
                    break;
                full_prob_line += " " + trim(next_line);
            }

            size_t start_paren = full_prob_line.find('(');
            size_t end_paren = full_prob_line.find(')');
            size_t pipe_pos = full_prob_line.find('|');

            if (start_paren == string::npos || end_paren == string::npos)
                continue;

            string prob_content = full_prob_line.substr(start_paren + 1, end_paren - start_paren - 1);
            stringstream prob_ss(prob_content);

            prob_ss >> node_name;

            int index = BayesNet.get_index(node_name);
            Graph_Node &node = BayesNet.get_node_by_index(index); // Get by reference

            vector<int> parent_indices; // Store indices, not names

            if (pipe_pos != string::npos && pipe_pos < end_paren)
            {
                string parents_str = full_prob_line.substr(pipe_pos + 1, end_paren - pipe_pos - 1);
                stringstream parent_ss(parents_str);
                string parent_name;

                while (parent_ss >> parent_name)
                {
                    if (parent_name.back() == ',')
                    {
                        parent_name = parent_name.substr(0, parent_name.length() - 1);
                    }
                    int parent_index = BayesNet.get_index(parent_name);
                    parent_indices.push_back(parent_index);

                    Graph_Node &parent_node = BayesNet.get_node_by_index(parent_index);
                    parent_node.add_child(index);
                }
            }

            node.set_Parent_Indices(parent_indices);

            vector<float> cpt;
            while (getline(myfile, line))
            {
                line = trim(line);
                if (line == "};")
                    break;
                if (line.empty())
                    continue;

                size_t close_paren = line.find(')');
                string prob_part;

                if (close_paren != string::npos)
                {
                    prob_part = line.substr(close_paren + 1);
                }
                else if (line.find("table") != string::npos)
                {
                    size_t table_pos = line.find("table");
                    prob_part = line.substr(table_pos + 5);
                }
                else
                {
                    prob_part = line;
                }

                stringstream ss_prob(prob_part);
                string token;
                while (ss_prob >> token)
                {
                    while (!token.empty() && (token.back() == ',' || token.back() == ';'))
                    {
                        token = token.substr(0, token.length() - 1);
                    }
                    if (!token.empty() && (isdigit(token[0]) || token[0] == '.' || token[0] == '-'))
                    {
                        cpt.push_back(atof(token.c_str()));
                    }
                }
            }
            node.set_CPT(cpt);
        }
    }

    myfile.close();
    return BayesNet;
}

// --- BIF Writer (Refactored) ---

void write_network(const char *filename, network &BayesNet)
{
    ofstream outfile(filename);

    if (!outfile.is_open())
    {
        cout << "Error: Could not open file " << filename << " for writing" << endl;
        return;
    }

    outfile << "// Bayesian Network (Learned by EM)" << endl
            << endl;

    int N = BayesNet.netSize();

    for (int i = 0; i < N; i++)
    {
        Graph_Node &node = BayesNet.get_node_by_index(i); // Get by ref

        outfile << "variable " << node.get_name() << " {" << endl;
        outfile << "  type discrete [ " << node.get_nvalues() << " ] = { ";

        vector<string> vals = node.get_values();
        for (int j = 0; j < (int)vals.size(); j++)
        {
            outfile << vals[j];
            if (j < (int)vals.size() - 1)
                outfile << ", ";
        }
        outfile << " };" << endl;
        outfile << "}" << endl;
    }

    // Use 4 decimal places as requested
    outfile << std::fixed << std::setprecision(4);
    for (int i = 0; i < N; i++)
    {
        Graph_Node &node = BayesNet.get_node_by_index(i);
        vector<int> parent_indices = node.get_Parent_Indices();
        vector<string> values = node.get_values();
        vector<float> cpt = node.get_CPT();

        outfile << "probability ( " << node.get_name();
        if (!parent_indices.empty())
        {
            outfile << " | ";
            for (int j = 0; j < (int)parent_indices.size(); j++)
            {
                outfile << BayesNet.get_node_by_index(parent_indices[j]).get_name();
                if (j < (int)parent_indices.size() - 1)
                    outfile << ", ";
            }
        }
        outfile << " ) {" << endl;

        vector<int> radices;
        radices.reserve(parent_indices.size());
        for (int p_index : parent_indices)
        {
            radices.push_back(BayesNet.get_node_by_index(p_index).get_nvalues());
        }

        int parent_combinations = 1;
        for (int r : radices)
            parent_combinations *= r;

        int cpt_index = 0;

        if (parent_indices.empty())
        {
            outfile << "    table ";
            for (int k = 0; k < (int)values.size(); k++)
            {
                if (cpt_index < (int)cpt.size())
                    outfile << cpt[cpt_index++];
                else
                    outfile << "-1"; // Should not happen after EM
                if (k < (int)values.size() - 1)
                    outfile << ", ";
            }
            outfile << ";" << endl;
        }
        else
        {
            for (int comb = 0; comb < parent_combinations; comb++)
            {
                vector<int> idx(parent_indices.size(), 0);
                int tmp = comb;
                for (int p = (int)parent_indices.size() - 1; p >= 0; p--)
                {
                    idx[p] = tmp % radices[p];
                    tmp /= radices[p];
                }

                outfile << "    ( ";
                for (int p = 0; p < (int)parent_indices.size(); p++)
                {
                    Graph_Node &pnode = BayesNet.get_node_by_index(parent_indices[p]);
                    auto pvals = pnode.get_values();
                    int vidx = idx[p];
                    outfile << pvals[vidx];
                    if (p < (int)parent_indices.size() - 1)
                        outfile << ", ";
                }
                outfile << " ) ";

                for (int k = 0; k < (int)values.size(); k++)
                {
                    if (cpt_index < (int)cpt.size())
                        outfile << cpt[cpt_index++];
                    else
                        outfile << "-1";
                    if (k < (int)values.size() - 1)
                        outfile << ", ";
                }
                outfile << ";" << endl;
            }
        }
        outfile << "};" << endl
                << endl;
    }

    outfile.close();
    // cout << "Network written to file: " << filename << endl; // Silenced for autograder
}

// --- New Functions for Data Reading and EM ---

/**
 * Reads the .dat file into a 2D vector of strings.
 * *** THIS IS THE FIXED VERSION ***
 */
vector<vector<string>> read_data(const char *filename)
{
    vector<vector<string>> records;
    ifstream file(filename);
    string line;

    if (!file.is_open())
    {
        cout << "Error: Could not open data file " << filename << endl;
        return records;
    }

    while (getline(file, line))
    {
        line = trim(line);
        if (line.empty())
            continue;

        vector<string> record;
        stringstream ss(line);
        string value;

        // --- THIS IS THE FIX ---
        // We are now splitting by commas (',') instead of spaces.
        while (getline(ss, value, ','))
        {

            // Remove quotes if they exist
            if (!value.empty() && value.front() == '"')
            {
                value = value.substr(1);
            }
            if (!value.empty() && value.back() == '"')
            {
                value = value.substr(0, value.length() - 1);
            }

            record.push_back(value);
        }
        // --- END OF FIX ---

        // Validate record size
        if (record.size() == 56)
        {
            records.push_back(record);
        }
        else if (!record.empty())
        {
            // This warning will now only trigger for truly bad lines
            cout << "Warning: Skipping malformed record with " << record.size() << " fields, expected 56." << endl;
        }
    }
    file.close();
    return records;
}

/**
 * Initializes all unknown (-1) probabilities to a uniform distribution.
 */
void initialize_probabilities(network &BayesNet)
{
    srand(time(0)); // Seed random number generator
    for (int i = 0; i < BayesNet.netSize(); ++i)
    {
        Graph_Node &node = BayesNet.get_node_by_index(i);
        vector<float> &cpt = node.get_CPT();
        int n_values = node.get_nvalues();

        if (cpt.empty())
        { // Handle nodes with no CPT data
            int parent_combinations = 1;
            for (int p_idx : node.get_Parent_Indices())
            {
                parent_combinations *= BayesNet.get_node_by_index(p_idx).get_nvalues();
            }
            cpt.resize(parent_combinations * n_values, -1.0f);
        }

        // Iterate in chunks of size n_values (one chunk per parent configuration)
        for (int j = 0; j < (int)cpt.size(); j += n_values)
        {
            if (cpt[j] == -1.0f)
            {
                // Unknown probabilities. Initialize randomly and normalize.
                float sum = 0.0f;
                vector<float> random_vals(n_values);
                for (int k = 0; k < n_values; ++k)
                {
                    // Use a small non-zero value to avoid log(0)
                    random_vals[k] = (float)rand() / RAND_MAX + 0.001f;
                    sum += random_vals[k];
                }
                for (int k = 0; k < n_values; ++k)
                {
                    cpt[j + k] = random_vals[k] / sum;
                }
            }
        }
    }
}

/**
 * Helper function to get the flat CPT index for a node given a full data record.
 * This will be very useful in your EM implementation.
 */
int get_cpt_index(const network &BayesNet, int node_index, const vector<string> &record)
{
    const Graph_Node &node = BayesNet.get_node_by_index(node_index);
    int n_values = node.get_nvalues();

    // 1. Get the index for the node's own value
    int value_index = node.get_value_index(record[node_index]);
    if (value_index == -1)
        return -1; // Should not happen in a complete record

    // 2. Get the index for the parent configuration
    vector<int> parent_indices = node.get_Parent_Indices();
    if (parent_indices.empty())
    {
        return value_index;
    }

    int parent_config_index = 0;
    int multiplier = 1;

    // Loop from the last parent to the first (this is how BIF CPTs are structured)
    for (int i = (int)parent_indices.size() - 1; i >= 0; --i)
    {
        int p_index = parent_indices[i];
        const Graph_Node &parent_node = BayesNet.get_node_by_index(p_index);

        string parent_value = record[p_index];
        int parent_value_index = parent_node.get_value_index(parent_value);
        if (parent_value_index == -1)
            return -1; // Parent value invalid

        parent_config_index += parent_value_index * multiplier;
        multiplier *= parent_node.get_nvalues();
    }

    // 3. Combine them
    int final_index = (parent_config_index * n_values) + value_index;
    return final_index;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//               YOUR IMPLEMENTATION GOES HERE
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/**
 * Calculates the joint probability of a single, complete record.
 * You will need this for the E-step.
 * Uses log-probabilities to prevent numerical underflow.
 */
float calculate_joint_probability(const network &BayesNet, const vector<string> &record)
{
    // Hint: Loop through all nodes, get their CPT index, and multiply their probabilities.
    // Using log probabilities (and summing them) is safer.
    // float joint_prob = 1.0;
    // ...
    // return joint_prob;

    double log_prob_sum = 0.0;
    for (int i = 0; i < BayesNet.netSize(); ++i)
    {
        int cpt_idx = get_cpt_index(BayesNet, i, record);
        if (cpt_idx == -1)
        {
            cerr << "Error calculating joint prob: invalid record." << endl;
            return 0.0;
        }

        float prob = BayesNet.get_node_by_index(i).get_CPT()[cpt_idx];
        if (prob == 0.0f)
        {
            return 0.0; // This configuration is impossible
        }
        log_prob_sum += log(prob);
    }

    return exp(log_prob_sum);
}

/**
 * Runs the Expectation-Maximization algorithm.
 */
void run_em(network &BayesNet, const vector<vector<string>> &records, int max_iterations, float convergence_threshold)
{

    cout << "Starting EM algorithm... (Stub)" << endl;
    cout << "You need to implement this function!" << endl;
    cout << "Max Iterations: " << max_iterations << ", Convergence: " << convergence_threshold << endl;

    for (int iter = 0; iter < max_iterations; ++iter)
    {

        // --- E-Step: Calculate Expected Sufficient Statistics (Counts) ---

        // 1. Initialize expected counts for this iteration.
        map<int, vector<float>> expected_counts;
        for (int i = 0; i < BayesNet.netSize(); ++i)
        {
            // Get the size of the CPT for node 'i'
            int cpt_size = BayesNet.get_node_by_index(i).get_CPT().size();

            // Create a vector of that size, initializing all values to 1.0f
            // (This is the Laplace "add-1" smoothing)
            expected_counts[i] = vector<float>(cpt_size, 1.0f);
        }

        // 2. Iterate over all records
        for (const auto &record : records)
        {

            // Find the missing_index (where record[i] == "?")
            int missing_index = -1; // Initialize to -1 (meaning "not found")
            for (int i = 0; i < (int)record.size(); ++i)
            {
                if (record[i] == "?")
                {
                    missing_index = i;
                    break; // Found it, stop searching
                }
            }

            if (missing_index == -1)
            {
                // --- Case 1: Complete Record ---
                // Add 1.0 to the observed counts for each node
                for (int i = 0; i < BayesNet.netSize(); ++i)
                {
                    int cpt_index = get_cpt_index(BayesNet, i, record);
                    // Check for a valid index
                    if (cpt_index != -1)
                    {
                        expected_counts[i][cpt_index] += 1.0f;
                    }
                }
            }
            else
            {
                // --- Case 2: Incomplete Record (one '?') ---

                // a) Calculate P(Missing | Evidence)
                Graph_Node &missing_node = BayesNet.get_node_by_index(missing_index);
                vector<string> possible_values = missing_node.get_values();
                vector<float> posterior(possible_values.size());

                float total_prob = 0.0f;

                // Loop through each possible value (e.g., "True", then "False")
                for (int i = 0; i < (int)possible_values.size(); ++i)
                {
                    // Create a copy of the record
                    vector<string> temp_record = record;
                    // Fill in the "?" with the value we're testing
                    temp_record[missing_index] = possible_values[i];

                    // Calculate the probability of the *entire* record with this value
                    float joint_prob = calculate_joint_probability(BayesNet, temp_record);

                    // Store this probability
                    posterior[i] = joint_prob;
                    total_prob += joint_prob;
                }

                // b) Normalize to get the posterior probability P(Missing | Evidence)
                if (total_prob > 0)
                {
                    // This turns the joint probabilities into a real posterior
                    // (e.g., [0.03, 0.07] -> [0.3, 0.7])
                    for (int i = 0; i < (int)posterior.size(); ++i)
                    {
                        posterior[i] /= total_prob;
                    }
                }
                else
                {
                    // Edge case: all configurations had 0 probability.
                    // Assign a uniform probability to avoid dividing by zero.
                    for (int i = 0; i < (int)posterior.size(); ++i)
                    {
                        posterior[i] = 1.0f / posterior.size();
                    }
                }

                // c) Add fractional counts based on the posterior
                for (int i = 0; i < (int)possible_values.size(); ++i)
                {
                    float weight = posterior[i];

                    // Optimization: If this value has 0 probability, skip it.
                    if (weight == 0.0f)
                    {
                        continue;
                    }

                    // Create the "completed" record for this specific value
                    vector<string> temp_record = record;
                    temp_record[missing_index] = possible_values[i];

                    // Loop through ALL nodes and add this record's
                    // weighted contribution to their counts.
                    for (int j = 0; j < BayesNet.netSize(); ++j)
                    {
                        // Find the CPT index for node 'j' using this temp_record
                        int cpt_index = get_cpt_index(BayesNet, j, temp_record);

                        // Important: Check for a valid index
                        if (cpt_index != -1)
                        {
                            // This is the core of the E-step:
                            // Add the fractional count (the weight)
                            expected_counts[j][cpt_index] += weight;
                        }
                    }
                }
            }
        } // End of loop over records

        // --- M-Step: Maximize Parameters (Update CPTs) ---

        float max_param_change = 0.0f;

        for (int i = 0; i < BayesNet.netSize(); ++i)
        {
            Graph_Node &node = BayesNet.get_node_by_index(i);
            vector<float> &cpt = node.get_CPT();
            vector<float> &counts = expected_counts[i];
            int n_values = node.get_nvalues();

            // Iterate in chunks (one per parent configuration)
            for (int j = 0; j < (int)cpt.size(); j += n_values)
            {

                float sum_of_counts = 0.0f;
                for (int k = 0; k < n_values; ++k)
                {
                    sum_of_counts += counts[j + k];
                }

                // b) Normalize to get new probabilities
                for (int k = 0; k < n_values; ++k)
                {
                    float new_prob;
                    if (sum_of_counts == 0.0f)
                    {
                        // This should not happen with add-1 smoothing
                        // But as a fallback, assign uniform probability
                        new_prob = 1.0f / n_values;
                    }
                    else
                    {
                        // Calculate the new probability
                        new_prob = counts[j + k] / sum_of_counts;
                    }

                    // Check how much the probability changed
                    float change = abs(new_prob - cpt[j + k]);
                    if (change > max_param_change)
                    {
                        max_param_change = change;
                    }

                    // Update the CPT with the new probability
                    cpt[j + k] = new_prob;
                }
            }
        }

        // --- Convergence Check ---
        cout << "Iteration " << iter + 1 << ": Max param change = " << max_param_change << endl;
        if (max_param_change < convergence_threshold)
        {
            cout << "Convergence reached!" << endl;
            break;
        }

        // --- Remove this break once you implement the loop ---
        // cout << "Stub: Breaking after 0 iterations." << endl;
        // break;
    }

    cout << "EM algorithm finished." << endl;
}

// --- Main Driver ---

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cout << "Usage: ./run.sh <bif_file> <data_file>" << endl;
        return 1;
    }

    string bif_filename = argv[1];
    string data_filename = argv[2];
    string output_filename = "solved_hailfinder.bif";

    // 1. Read network structure
    network BayesNet = read_network(bif_filename.c_str());
    if (BayesNet.netSize() == 0)
    {
        cout << "Error: Failed to read network." << endl;
        return 1;
    }

    // 2. Read data
    vector<vector<string>> records = read_data(data_filename.c_str());
    if (records.empty())
    {
        cout << "Error: Failed to read data or data file is empty." << endl;
        return 1;
    }

    // 3. Initialize unknown probabilities
    initialize_probabilities(BayesNet);

    // 4. Run EM algorithm
    // 50 iterations or 1e-4 change is a reasonable default.
    run_em(BayesNet, records, 50, 1e-4);

    // 5. Write the solved network
    write_network(output_filename.c_str(), BayesNet);

    return 0;
}
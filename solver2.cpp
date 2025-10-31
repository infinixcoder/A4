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

class Graph_Node
{
private:
    string Node_Name;
    vector<int> Children;
    vector<int> Parent_Indices;
    int nvalues;
    vector<string> values;
    map<string, int> value_to_index;
    vector<float> CPT;

public:
    Graph_Node(string name, int n, vector<string> vals)
    {
        Node_Name = name;
        nvalues = n;
        values = vals;

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
    {
        return Parent_Indices;
    }

    vector<float> &get_CPT()
    {
        return CPT;
    }

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
            return -1;
        }
    }

    void set_CPT(vector<float> new_CPT)
    {
        CPT.clear();
        CPT = new_CPT;
    }

    void set_Parent_Indices(vector<int> Parent_Node_Indices)
    {
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
    vector<Graph_Node> nodes;
    map<string, int> name_to_index;

public:
    int addNode(Graph_Node node)
    {
        int index = nodes.size();
        nodes.push_back(node);
        name_to_index[node.get_name()] = index;
        return 0;
    }

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
            Graph_Node &node = BayesNet.get_node_by_index(index);

            vector<int> parent_indices;

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

void write_network(const char *filename, network &BayesNet)
{
    ofstream outfile(filename);

    if (!outfile.is_open())
    {
        cout << "Error: Could not open file " << filename << " for writing" << endl;
        return;
    }

    outfile << "" << endl
            << endl;

    int N = BayesNet.netSize();

    for (int i = 0; i < N; i++)
    {
        Graph_Node &node = BayesNet.get_node_by_index(i);

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
                    outfile << "-1";
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
}

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

            record.push_back(value);
        }

        if (record.size() == 56)
        {
            records.push_back(record);
        }
        else if (!record.empty())
        {

            cout << "Warning: Skipping malformed record with " << record.size() << " fields, expected 56." << endl;
        }
    }
    file.close();
    return records;
}

void initialize_probabilities(network &BayesNet)
{
    srand(time(0));
    for (int i = 0; i < BayesNet.netSize(); ++i)
    {
        Graph_Node &node = BayesNet.get_node_by_index(i);
        vector<float> &cpt = node.get_CPT();
        int n_values = node.get_nvalues();

        if (cpt.empty())
        {
            int parent_combinations = 1;
            for (int p_idx : node.get_Parent_Indices())
            {
                parent_combinations *= BayesNet.get_node_by_index(p_idx).get_nvalues();
            }
            cpt.resize(parent_combinations * n_values, -1.0f);
        }

        for (int j = 0; j < (int)cpt.size(); j += n_values)
        {
            if (cpt[j] == -1.0f)
            {

                float sum = 0.0f;
                vector<float> random_vals(n_values);
                for (int k = 0; k < n_values; ++k)
                {

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

int get_cpt_index(const network &BayesNet, int node_index, const vector<string> &record)
{
    const Graph_Node &node = BayesNet.get_node_by_index(node_index);
    int n_values = node.get_nvalues();

    int value_index = node.get_value_index(record[node_index]);
    if (value_index == -1)
        return -1;

    vector<int> parent_indices = node.get_Parent_Indices();
    if (parent_indices.empty())
    {
        return value_index;
    }

    int parent_config_index = 0;
    int multiplier = 1;

    for (int i = (int)parent_indices.size() - 1; i >= 0; --i)
    {
        int p_index = parent_indices[i];
        const Graph_Node &parent_node = BayesNet.get_node_by_index(p_index);

        string parent_value = record[p_index];
        int parent_value_index = parent_node.get_value_index(parent_value);
        if (parent_value_index == -1)
            return -1;

        parent_config_index += parent_value_index * multiplier;
        multiplier *= parent_node.get_nvalues();
    }

    int final_index = (parent_config_index * n_values) + value_index;
    return final_index;
}

float calculate_joint_probability(const network &BayesNet, const vector<string> &record)
{

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
            return 0.0;
        }
        log_prob_sum += log(prob);
    }

    return exp(log_prob_sum);
}

void run_em(network &BayesNet, const vector<vector<string>> &records, int max_iterations, float convergence_threshold)
{

    cout << "Starting EM algorithm... (Stub)" << endl;
    cout << "You need to implement this function!" << endl;
    cout << "Max Iterations: " << max_iterations << ", Convergence: " << convergence_threshold << endl;

    for (int iter = 0; iter < max_iterations; ++iter)
    {

        map<int, vector<float>> expected_counts;
        for (int i = 0; i < BayesNet.netSize(); ++i)
        {

            int cpt_size = BayesNet.get_node_by_index(i).get_CPT().size();

            expected_counts[i] = vector<float>(cpt_size, 0.6f);
        }

        for (const auto &record : records)
        {

            int missing_index = -1;
            for (int i = 0; i < (int)record.size(); ++i)
            {
                if (record[i] == "?")
                {
                    missing_index = i;
                    break;
                }
            }

            if (missing_index == -1)
            {

                for (int i = 0; i < BayesNet.netSize(); ++i)
                {
                    int cpt_index = get_cpt_index(BayesNet, i, record);

                    if (cpt_index != -1)
                    {
                        expected_counts[i][cpt_index] += 1.0f;
                    }
                }
            }
            else
            {

                Graph_Node &missing_node = BayesNet.get_node_by_index(missing_index);
                vector<string> possible_values = missing_node.get_values();
                vector<float> posterior(possible_values.size());

                float total_prob = 0.0f;

                for (int i = 0; i < (int)possible_values.size(); ++i)
                {

                    vector<string> temp_record = record;

                    temp_record[missing_index] = possible_values[i];

                    float joint_prob = calculate_joint_probability(BayesNet, temp_record);

                    posterior[i] = joint_prob;
                    total_prob += joint_prob;
                }

                if (total_prob > 0)
                {

                    for (int i = 0; i < (int)posterior.size(); ++i)
                    {
                        posterior[i] /= total_prob;
                    }
                }
                else
                {

                    for (int i = 0; i < (int)posterior.size(); ++i)
                    {
                        posterior[i] = 1.0f / posterior.size();
                    }
                }

                for (int i = 0; i < (int)possible_values.size(); ++i)
                {
                    float weight = posterior[i];

                    if (weight == 0.0f)
                    {
                        continue;
                    }

                    vector<string> temp_record = record;
                    temp_record[missing_index] = possible_values[i];

                    for (int j = 0; j < BayesNet.netSize(); ++j)
                    {

                        int cpt_index = get_cpt_index(BayesNet, j, temp_record);

                        if (cpt_index != -1)
                        {

                            expected_counts[j][cpt_index] += weight;
                        }
                    }
                }
            }
        }

        float max_param_change = 0.0f;

        for (int i = 0; i < BayesNet.netSize(); ++i)
        {
            Graph_Node &node = BayesNet.get_node_by_index(i);
            vector<float> &cpt = node.get_CPT();
            vector<float> &counts = expected_counts[i];
            int n_values = node.get_nvalues();

            for (int j = 0; j < (int)cpt.size(); j += n_values)
            {

                float sum_of_counts = 0.0f;
                for (int k = 0; k < n_values; ++k)
                {
                    sum_of_counts += counts[j + k];
                }

                for (int k = 0; k < n_values; ++k)
                {
                    float new_prob;
                    if (sum_of_counts == 0.0f)
                    {

                        new_prob = 1.0f / n_values;
                    }
                    else
                    {

                        new_prob = counts[j + k] / sum_of_counts;
                    }

                    float change = abs(new_prob - cpt[j + k]);
                    if (change > max_param_change)
                    {
                        max_param_change = change;
                    }

                    cpt[j + k] = new_prob;
                }
            }
        }

        cout << "Iteration " << iter + 1 << ": Max param change = " << max_param_change << endl;
        if (max_param_change < convergence_threshold)
        {
            cout << "Convergence reached!" << endl;
            break;
        }
    }

    cout << "EM algorithm finished." << endl;
}

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

    network BayesNet = read_network(bif_filename.c_str());
    if (BayesNet.netSize() == 0)
    {
        cout << "Error: Failed to read network." << endl;
        return 1;
    }

    vector<vector<string>> records = read_data(data_filename.c_str());
    if (records.empty())
    {
        cout << "Error: Failed to read data or data file is empty." << endl;
        return 1;
    }

    initialize_probabilities(BayesNet);

    run_em(BayesNet, records, 50, 2e-7);

    write_network(output_filename.c_str(), BayesNet);

    return 0;
}
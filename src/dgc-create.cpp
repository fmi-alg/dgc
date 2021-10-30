#include <iostream>
#include <vector>
#include <functional>
#include <boost/serialization/utility.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include "c.h"
#include "Timer.h"
#include "ch.h"
#include <fstream>
#include <omp.h>
#include <getopt.h>
#include <chrono>
#include <ctime>
#include <sstream>

using namespace std;

typedef pair<Rank, Distance> LabelElement;

string getNameOfQueryType(int query_type)
{
    string query_name;
    if (query_type == c::QUERY_BUCKETS)
    {
        query_name = "buckets";
    }
    else if (query_type == c::QUERY_CH)
    {
        query_name = "ch";
    }
    else if (query_type == c::QUERY_PQ)
    {
        query_name = "pq";
    }
    else if (query_type == c::QUERY_CHHL)
    {
        query_name = "chhl";
    }
    else if (query_type == c::QUERY_HL)
    {
        query_name = "hl";
    }
    else if (query_type == c::QUERY_BUCKETSMAPS)
    {
        query_name = "bucketsmaps";
    }
    else if (query_type == c::QUERY_PQSoD)
    {
        query_name = "PQ with stall on demand";
    }
    else
    {
        query_name = "unknown";
    }
    return query_name;
}

string getTime()
{
    auto time1 = std::chrono::system_clock::now();
    time_t time2 = chrono::system_clock::to_time_t(time1);
    return ctime(&time2);
}

struct Measurement
{
    double time;
    int query_type;
    string file_name;
    string today;
    int num_runs;
    double average_label_size;
    double average_label_size_fwd;
    long num_nodes;
    long total_num_hubs;
    vector<long> space_per_group;
    vector<Rank> group_bounds;
    vector<int> query_analysis;
    string get_output()
    {
        string output = "";
        output += "Input file: " + file_name + "\n";
        output += "Date: " + today;
        output += "Query type: " + getNameOfQueryType(query_type) + "\n";
        output += "Num runs: " + to_string(num_runs) + "\n";
        output += "Time in millis: " + to_string(time * 1000) + "\n";
        output += "Average label size: " + to_string(average_label_size) + "\n";
        output += "Average label size fwd: " + to_string(average_label_size_fwd) + "\n";
        output += "Total number of edges: " + to_string(total_num_hubs) + "\n";
        output += "Total number of nodes: " + to_string(num_nodes) + "\n";
        output += "Bounds %:";
        for (Group i = 0; i < group_bounds.size() - 1; i++)
        {
            output += " " + to_string((double)(group_bounds[i]) / group_bounds.at(group_bounds.size() - 1));
        }
        output += "\n";
        output += "Size per label: " + to_string(sizeof(LabelElement)) + " Bytes\n";
        output += "Size per rank: " + to_string(sizeof(Rank)) + " Bytes\n";
        double total_space_disk = total_num_hubs * sizeof(LabelElement) + 2 * sizeof(Rank) * num_nodes;
        double total_space_ram = total_space_disk + 2 * sizeof(Distance) * num_nodes;
        output += "Total Space Disk: " + to_string(total_space_disk / (1024 * 1024)) + " MB, " + to_string(total_space_disk / (1024 * 1024 * 1024)) + " GB" + "\n";
        output += "Total Space RAM: " + to_string(total_space_ram / (1024 * 1024)) + " MB, " + to_string(total_space_ram / (1024 * 1024 * 1024)) + " GB" + "\n";
        output += "Space per group:";
        for (long s : space_per_group)
        {
            output += " " + to_string((double)(s) / total_num_hubs);
        }
        output += "\n";
        for (int i = 0; i < query_analysis.size(); i++)
        {
            if (i == 0)
            {
                output += "Average number of labels per query: " + to_string(query_analysis[i]) + "\n";
            }
            else if (i == 1)
            {
                output += "Average number of edges per query: " + to_string(query_analysis[i]) + "\n";
            }
        }
        return output;
    }
};

class KHopLabeling
{
public:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &_fwdLabels;
        ar &_bwdLabels;
        ar &_groupBounds;
    }
    KHopLabeling()
    {
    }

    void constructFromCH(string ch_binary, vector<Level> group_bounds_in_levels, int num_threads = c::NO_ENTRY, int prune_labels = c::NO_ENTRY)
    {
        _numHops = group_bounds_in_levels.size() + 1;
        Timer timer;
        timer.start();
        _constructLabels(ch_binary, group_bounds_in_levels, num_threads, prune_labels);
        timer.stop();
        cout << "Finished construction in " << timer.secs() << " seconds" << endl;
        _fwd_distances_per_thread.clear();
        _bwd_distances_per_thread.clear();
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
    }

    void constructHL(string file, int nr_of_threads = c::NO_ENTRY)
    {
        _numHops = 1;
        initializeWithCH(file);
        _numNodes = _ch.nofNodes();
        _fwdLabels.resize(_numNodes);
        _bwdLabels.resize(_numNodes);
        if (nr_of_threads == c::NO_ENTRY)
        {
            nr_of_threads = omp_get_max_threads();
        }
        cout << "Start constructing labels with " << nr_of_threads << " threads." << endl;
        Timer construction_timer;
        construction_timer.start();
        _fwd_distances_per_thread.resize(nr_of_threads);
        _bwd_distances_per_thread.resize(nr_of_threads);
        for (int i = 0; i < nr_of_threads; i++)
        {
            _fwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
            _bwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
        }
        _groupBounds.clear();
        _groupBounds.push_back(_numNodes);
        {
            Rank i = 0;
            int round = 0;
            long next_print_step = _numNodes / 30;
            long next_print = next_print_step;
            while (i < _numNodes)
            {
                int first_node = i;
                i = _numNodes - _ch.getFirstNodeOfLevel(_ch.getMaxLevel() - round);
                round++;
#pragma omp parallel for num_threads(nr_of_threads)
                for (Rank x = first_node; x < i; x++)
                {
                    int tid = omp_get_thread_num();
                    Rank rank = _numNodes - 1 - x;
                    _constructLabelHL(rank, tid);
                }
            }
        }
        construction_timer.stop();
        cout << "Finished construciton in " << construction_timer.secs() << " seconds." << endl;
        _fwd_distances_per_thread.clear();
        _bwd_distances_per_thread.clear();
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
    }

    void constructCHHL(string file, Level level_bound, int nr_of_threads = c::NO_ENTRY)
    {
        Timer construction_timer;
        construction_timer.start();
        _numHops = 2;
        vector<Rank> first_node_of_level;
        Level max_level;
        {
            CH ch;
            ch.readFromBinary(file);
            max_level = ch.getMaxLevel();
            _numNodes = ch.nofNodes();
            _fwdLabels.resize(_numNodes);
            _bwdLabels.resize(_numNodes);
            first_node_of_level.resize(ch.getMaxLevel() + 1);
            for (int i = 0; i <= ch.getMaxLevel(); i++)
            {
                first_node_of_level.at(i) = ch.getFirstNodeOfLevel(i);
            }
            cout << "Working with " << nr_of_threads << " threads. Start construction..." << endl;
#pragma omp parallel for num_threads(nr_of_threads)
            for (Rank rank = 0; rank < _numNodes; rank++)
            {
                _fwdLabels.at(rank) = ch.getNeighbors(rank, true);
                _bwdLabels.at(rank) = ch.getNeighbors(rank, false);
            }
            cout << "low ranks finished" << endl;
        }

        Rank node_bound = first_node_of_level[level_bound];
        if (nr_of_threads == c::NO_ENTRY)
        {
            nr_of_threads = 1;
        }

        _fwd_distances_per_thread.resize(nr_of_threads);
        _bwd_distances_per_thread.resize(nr_of_threads);
        for (int i = 0; i < nr_of_threads; i++)
        {
            _fwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
            _bwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
        }
        _groupBounds.clear();
        _groupBounds.push_back(node_bound);
        _groupBounds.push_back(_numNodes);
        {
            Rank i = 0;
            int round = 0;
            while (i < _numNodes - node_bound)
            {
                int first_node = i;
                round++;
                i = _numNodes - first_node_of_level[max_level - round];

#pragma omp parallel for num_threads(nr_of_threads)
                for (Rank x = first_node; x < i; x++)
                {
                    int tid = omp_get_thread_num();
                    Rank rank = _numNodes - 1 - x;
                    _constructLabelWithoutCH(rank, tid);
                }
            }
        }
        construction_timer.stop();
        for (int i = 0; i < nr_of_threads; i++)
        {
            _fwd_distances_per_thread[i].clear();
            _bwd_distances_per_thread[i].clear();
        }
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        cout << "Finished. Construciton took " << construction_timer.secs() << " seconds." << endl;
    }

    long getSizeOfLabeling()
    {
        long total_num_hubs = 0;
        for (long i = 0; i < _numNodes; i++)
        {
            total_num_hubs += _fwdLabels[i].size() + _bwdLabels[i].size();
        }
        return total_num_hubs * sizeof(LabelElement) + 2 * sizeof(Rank) * _numNodes;
    }

    vector<Rank> getBounds()
    {
        return _groupBounds;
    }

    vector<Level> sampleBestBoundsWithEstimation(string ch_binary, int num_hops, int num_threads = c::NO_ENTRY)
    {
        cout << "Start sampling with " << (num_threads == c::NO_ENTRY ? 1 : num_threads) << " thread(s)" << endl;
        cout << "Time: " << getTime();
        Timer timer;
        timer.start();
        vector<Level> bounds;
        vector<Level> best_product_solution;
        double best_time = 999999999;
        double best_space = 99999999999;
        double best_product;
        vector<Level> max_levels;
        {
            Level tmp = 25 / (num_hops - 1);
            for (int i = 0; i < num_hops - 2; i++)
            {
                max_levels.push_back(tmp);
                tmp += 25 / (num_hops - 1);
            }
            max_levels.push_back(40);
        }

        vector<bool> vacant;
        {
            long tmp = 1;
            for (int i = 0; i < num_hops - 1; i++)
            {
                tmp *= max_levels.at(i);
            }
            vacant.assign(tmp, true);
        }

        for (int i = 0; i < num_hops - 1; i++)
        {
            int level = rand() % max_levels.at(i);
            bounds.push_back(level + 1);
        }
        sort(bounds.begin(), bounds.end());
        {
            Rank index = 0;
            long p = 1;
            for (int i = 0; i < num_hops - 1; i++)
            {
                index += (bounds.at(i) - 1) * p;
                p *= max_levels.at(i);
            }
            vacant.at(index) = false;
        }
        initializeWithCH(ch_binary);
        cout << "Initial bounds:";
        for (int i = 0; i < bounds.size(); i++)
        {
            cout << " " << bounds.at(i);
        }
        cout << endl;
        int num_runs = 100;
        vector<double> result = estimateQueryComplexityWithCH(bounds, num_runs, num_threads);
        //DEBUG
        cout << "result:" << endl;
        cout << "edges: " << result[0] << endl;
        cout << "nodes: " << result[1] << endl;
        cout << "labels: " << result[2] << endl;
        cout << "space: " << result[3] << endl;
        cout << "time: " << result[0] + result[1] << endl;
        cout << "product: " << (result[0] + result[1]) * result[3] << endl;
        //DEBUG END
        best_time = result[0] + result[1];
        best_space = result[3];
        best_product = best_time * best_space;
        best_product_solution = bounds;
        int num_queries = 1000;
        if (num_queries > vacant.size())
        {
            num_queries = vacant.size();
        }
        int query_counter = 1;
        for (int q = 1; q < num_queries; q++)
        {
            bool is_random = (rand() % num_queries) > q;
            bool found = false;
            int counter = 0;
            int max_rounds = 10000;
            if (is_random)
            {
                cout << "Preparing random query." << endl;
                while (!found && counter < max_rounds)
                {
                    counter++;
                    for (int i = 0; i < num_hops - 1; i++)
                    {
                        Level level = 1 + (rand() % (max_levels.at(i) - 1));
                        bounds.at(i) = level;
                    }
                    sort(bounds.begin(), bounds.end());
                    {
                        Rank index = 0;
                        long p = 1;
                        for (int i = 0; i < num_hops - 1; i++)
                        {
                            index += (bounds.at(i) - 1) * p;
                            p *= max_levels.at(i);
                        }
                        found = vacant.at(index);
                        for (int i = 0; i < bounds.size() - 1; i++)
                        {
                            if (bounds.at(i) == bounds.at(i + 1))
                            {
                                found = false;
                            }
                        }
                    }
                }
            }
            else
            {
                cout << "Prepare query close to optimal solution for product." << endl;
                bounds = best_product_solution;
                Rank index = 0;
                {
                    long p = 1;
                    for (int i = 0; i < num_hops - 1; i++)
                    {
                        index += (bounds.at(i) - 1) * p;
                        p *= max_levels.at(i);
                    }
                }
                while (!found && counter < max_rounds)
                {
                    counter++;
                    int bound_index = rand() % bounds.size();
                    int dif = 2 * (rand() % 2) - 1;
                    bounds.at(bound_index) += dif;
                    if (bounds.at(bound_index) <= 0)
                    {
                        bounds.at(bound_index) = 1;
                    }
                    if (bounds.at(bound_index) >= max_levels.at(bound_index))
                    {
                        bounds.at(bound_index) = max_levels.at(bound_index) - 1;
                    }
                    if (bound_index > 0 && bounds.at(bound_index) <= bounds.at(bound_index - 1))
                    {
                        bounds.at(bound_index) = bounds.at(bound_index - 1) + 1;
                    }
                    if (bound_index < bounds.size() - 1 && bounds.at(bound_index) >= bounds.at(bound_index + 1))
                    {
                        bounds.at(bound_index) = bounds.at(bound_index + 1) - 1;
                    }
                    index = 0;
                    {
                        long p = 1;
                        for (int i = 0; i < num_hops - 1; i++)
                        {
                            index += (bounds.at(i) - 1) * p;
                            p *= max_levels.at(i);
                        }
                    }
                    found = vacant.at(index);
                }
            }
            if (!found)
            {
                break;
            }
            //DEBUG
            cout << "Checking bounds";
            for (int i = 0; i < bounds.size(); i++)
            {
                cout << " " << bounds[i];
            }
            cout << endl;
            //DEBUG END
            {
                Rank index = 0;
                long p = 1;
                for (int i = 0; i < num_hops - 1; i++)
                {
                    index += (bounds.at(i) - 1) * p;
                    p *= max_levels.at(i);
                }
                assert(vacant.at(index));
                vacant.at(index) = false;
            }
            result = estimateQueryComplexityWithCH(bounds, num_runs, num_threads);
            query_counter++;

            double this_space = result[3];
            double this_time = result[0] + result[1];
            if (this_time * this_space < best_product)
            {
                best_product = this_time * this_space;
                best_product_solution = bounds;
            }
            cout << "num queries: " << query_counter << endl;

            //DEBUG
            cout << "result:" << endl;
            cout << "edges: " << result[0] << endl;
            cout << "nodes: " << result[1] << endl;
            cout << "labels: " << result[2] << endl;
            cout << "space: " << result[3] << endl;
            cout << "time: " << result[0] + result[1] << endl;
            cout << "product: " << (result[0] + result[1]) * result[3] << endl;
            //DEBUG END
            cout << endl;
            cout << "best results so far:" << endl;
            cout << "best product: " << best_product;
            cout << " with bounds";
            for (int i = 0; i < best_product_solution.size(); i++)
            {
                cout << " " << best_product_solution.at(i);
            }
            cout << endl;
            cout << "progress: " << q * 100 / num_queries << "%" << endl;
            cout << endl;
            //DEBUG END
        }
        timer.stop();
        //DEBUG
        cout << endl;
        cout << "Finished." << endl;
        cout << "Num queries: " << query_counter << endl;
        cout << "Runs per query: " << num_runs << endl;
        cout << "Total time in seconds: " << timer.secs() << endl;
        cout << "Time per run in seconds: " << timer.secs() / (query_counter * num_runs) << endl;
        //DEBUG END
        return best_product_solution;
    }

    void spaceStressTest()
    {
        vector<vector<Rank>> all_vectors;
        while (true)
        {
            vector<Rank> new_vector(_numNodes);
            all_vectors.push_back(new_vector);
            cout << "Size of vector: " << all_vectors.size() << endl;
        }
    }

    vector<double> estimateQueryComplexityWithCH(vector<Level> level_bounds, int num_runs = 100, int nr_of_threads = c::NO_ENTRY)
    {
        _groupBounds.clear();
        for (int i = 0; i < level_bounds.size(); i++)
        {
            _groupBounds.push_back(_ch.getFirstNodeOfLevel(level_bounds[i]));
        }
        _groupBounds.push_back(_numNodes);
        _numHops = _groupBounds.size();

        if (nr_of_threads == c::NO_ENTRY)
        {
            nr_of_threads = 1;
        }

        vector<vector<long>> sum_edges_per_group_and_thread(_groupBounds.size());
        vector<vector<long>> sum_nodes_per_group_and_thread(_groupBounds.size());
        vector<vector<long>> sum_labels_per_group_and_thread(_groupBounds.size());
        vector<vector<long>> sum_label_size_per_group_and_thread(_groupBounds.size());

        _ch.setNumberOfThreads(nr_of_threads);

        if (_fwd_distances_per_thread.size() != nr_of_threads)
        {
            _fwd_distances_per_thread.resize(nr_of_threads);
            for (int i = 0; i < nr_of_threads; i++)
            {
                _fwd_distances_per_thread.at(i).resize(_numNodes, c::NO_ENTRY);
            }
        }

        for (int k = 0; k < _groupBounds.size(); k++)
        {
            Rank first_rank = k == 0 ? 0 : _groupBounds[k - 1] + 1;
            Rank last_rank = _groupBounds[k];
            sum_edges_per_group_and_thread[k].assign(nr_of_threads, 0);
            sum_nodes_per_group_and_thread[k].assign(nr_of_threads, 0);
            sum_labels_per_group_and_thread[k].assign(nr_of_threads, 0);
            sum_label_size_per_group_and_thread[k].assign(nr_of_threads, 0);
#pragma omp parallel for num_threads(nr_of_threads)
            for (int i = 0; i < num_runs; i++)
            {
                int t_id = omp_get_thread_num();
                Rank source = first_rank + (rand() % (last_rank - first_rank));
                vector<long> res = getQueryComplexityWithCH(source, t_id);
                sum_edges_per_group_and_thread[k][t_id] += res[0];
                sum_nodes_per_group_and_thread[k][t_id] += res[1];
                sum_labels_per_group_and_thread[k][t_id] += res[2];
                sum_label_size_per_group_and_thread[k][t_id] += res[3];
            }
        }

        double sum_edges = 0;
        double sum_nodes = 0;
        double sum_labels = 0;
        double sum_label_size = 0;
        for (int k = 0; k < _groupBounds.size(); k++)
        {
            Rank first_rank = k == 0 ? 0 : _groupBounds[k - 1] + 1;
            Rank last_rank = _groupBounds[k];
            double ratio = (double)(last_rank - first_rank + 1) / _numNodes;
            long sum_edges_per_group = 0;
            long sum_nodes_per_group = 0;
            long sum_labels_per_group = 0;
            long sum_label_size_per_group = 0;
            for (int i = 0; i < nr_of_threads; i++)
            {
                sum_edges_per_group += sum_edges_per_group_and_thread[k][i];
                sum_nodes_per_group += sum_nodes_per_group_and_thread[k][i];
                sum_labels_per_group += sum_labels_per_group_and_thread[k][i];
                sum_label_size_per_group += sum_label_size_per_group_and_thread[k][i];
            }
            sum_edges += sum_edges_per_group * ratio;
            sum_nodes += sum_nodes_per_group * ratio;
            sum_labels += sum_labels_per_group * ratio;
            sum_label_size += sum_label_size_per_group * ratio;
        }
        vector<double> result;
        result.push_back(sum_edges / num_runs);
        result.push_back(sum_nodes / num_runs);
        result.push_back(sum_labels / num_runs);
        result.push_back(sum_label_size / num_runs);
        return result;
    }

    void readFromBinary(string file)
    {
        cout << "Start binary read" << endl;
        ifstream myHLfile(file, ios_base::in | ios_base::binary);
        boost::archive::binary_iarchive ia(myHLfile);
        ia &(*this);
        _numNodes = _fwdLabels.size();
        _numHops = _groupBounds.size();
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        _added_by_fwd.assign(_numNodes, c::NO_ENTRY);
        _added_by_bwd.assign(_numNodes, c::NO_ENTRY);
        cout << "Finished reading." << endl;
    }

    long getNrOfNodes()
    {
        return _numNodes;
    }

    long getNrOfFwdEdges()
    {
        long fwd_edges = 0;
        for (long i = 0; i < _numNodes; i++)
        {
            fwd_edges += _fwdLabels[i].size();
        }
        return fwd_edges;
    }

    long getNrOfBwdEdges()
    {
        long bwd_edges = 0;
        for (long i = 0; i < _numNodes; i++)
        {
            bwd_edges += _bwdLabels[i].size();
        }
        return bwd_edges;
    }

    Distance getDistanceCHHL(Rank source, Rank target, int thread_id = c::NO_ENTRY)
    {
        if (source == target)
        {
            return 0;
        }
        vector<Distance> &fwdDistances = thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread.at(thread_id);
        vector<Distance> &bwdDistances = thread_id == c::NO_ENTRY ? _bwdDistances : _bwd_distances_per_thread.at(thread_id);
        priority_queue<pair<Distance, Rank>, vector<pair<Distance, Rank>>, std::greater<pair<Distance, Rank>>> pQ;
        vector<Rank> fwd_visited_nodes;
        vector<Rank> bwd_visited_nodes;
        pQ.push(make_pair(0, source));
        fwdDistances[source] = 0;
        fwd_visited_nodes.push_back(source);
        while (!pQ.empty())
        {
            pair<Distance, Rank> next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (fwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (LabelElement &label : _bwdLabels[rank])
                {
                    if (fwdDistances[label.first] != c::NO_ENTRY)
                    {
                        Distance new_dist = fwdDistances[label.first] + label.second;
                        if (new_dist < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    Group group = _getGroupOfNode(rank);
                    for (LabelElement &label : _fwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (fwdDistances[label.first] == c::NO_ENTRY)
                        {
                            fwdDistances[label.first] = new_distance;
                            fwd_visited_nodes.push_back(label.first);
                            if (group == 0)
                            {
                                pQ.push(make_pair(new_distance, label.first));
                            }
                        }
                        else if (new_distance < fwdDistances[label.first])
                        {
                            fwdDistances[label.first] = new_distance;
                            if (group == 0)
                            {
                                pQ.push(make_pair(new_distance, label.first));
                            }
                        }
                    }
                }
            }
        }
        pQ.push(make_pair(0, target));
        bwdDistances[target] = 0;
        bwd_visited_nodes.push_back(target);
        while (!pQ.empty())
        {
            pair<Distance, Rank> next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (bwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (LabelElement &label : _fwdLabels[rank])
                {
                    if (bwdDistances[label.first] != c::NO_ENTRY)
                    {
                        Distance new_dist = bwdDistances[label.first] + label.second;
                        if (new_dist < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    Group group = _getGroupOfNode(rank);
                    for (LabelElement &label : _bwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (bwdDistances[label.first] == c::NO_ENTRY)
                        {
                            bwdDistances[label.first] = new_distance;
                            bwd_visited_nodes.push_back(label.first);
                            if (group == 0)
                            {
                                pQ.push(make_pair(new_distance, label.first));
                            }
                        }
                        else if (new_distance < bwdDistances[label.first])
                        {
                            bwdDistances[label.first] = new_distance;
                            if (group == 0)
                            {
                                pQ.push(make_pair(new_distance, label.first));
                            }
                        }
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        for (Rank rank : fwd_visited_nodes)
        {
            if (bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = bwdDistances[rank] + fwdDistances[rank];
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                }
            }
            fwdDistances[rank] = c::NO_ENTRY;
        }
        for (Rank rank : bwd_visited_nodes)
        {
            bwdDistances[rank] = c::NO_ENTRY;
        }

        return best_dist;
    }

    Distance getDistanceBuckets(Rank source, Rank target)
    {
        if (source == target)
        {
            return 0;
        }
        vector<vector<Rank>> fwd_buckets(_numHops);
        vector<vector<Rank>> bwd_buckets(_numHops);
        vector<Rank> fwd_visited_nodes;
        vector<Rank> bwd_visited_nodes;
        fwd_visited_nodes.reserve(3000);
        bwd_visited_nodes.reserve(3000);
        Group source_group = _getGroupOfNode(source);
        _fwdDistances[source] = 0;
        fwd_visited_nodes.push_back(source);
        for (LabelElement &hub : _fwdLabels[source])
        {
            Distance new_distance = hub.second;
            if (_fwdDistances[hub.first] == c::NO_ENTRY)
            {
                _fwdDistances[hub.first] = new_distance;
                fwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (source_group < next_group)
                {
                    fwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < _fwdDistances[hub.first])
            {
                _fwdDistances[hub.first] = new_distance;
            }
        }
        for (Group group = source_group + 1; group < _numHops; group++)
        {
            for (Rank rank : fwd_buckets[group])
            {
                Distance distance = _fwdDistances[rank];
                for (LabelElement &hub : _fwdLabels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (_fwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        _fwdDistances[hub.first] = new_distance;
                        fwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            fwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < _fwdDistances[hub.first])
                    {
                        _fwdDistances[hub.first] = new_distance;
                    }
                }
            }
        }
        Group target_group = _getGroupOfNode(target);
        _bwdDistances[target] = 0;
        bwd_visited_nodes.push_back(target);
        for (LabelElement &hub : _bwdLabels[target])
        {
            Distance new_distance = hub.second;
            if (_bwdDistances[hub.first] == c::NO_ENTRY)
            {
                _bwdDistances[hub.first] = new_distance;
                bwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (target_group < next_group)
                {
                    bwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < _bwdDistances[hub.first])
            {
                _bwdDistances[hub.first] = new_distance;
            }
        }
        for (Group group = target_group + 1; group < _numHops; group++)
        {
            for (Rank rank : bwd_buckets[group])
            {
                Distance distance = _bwdDistances[rank];
                for (LabelElement &hub : _bwdLabels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (_bwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        _bwdDistances[hub.first] = new_distance;
                        bwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            bwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < _bwdDistances[hub.first])
                    {
                        _bwdDistances[hub.first] = new_distance;
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        for (Rank rank : fwd_visited_nodes)
        {
            if (_bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = _bwdDistances[rank] + _fwdDistances[rank];
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                }
            }
            _fwdDistances[rank] = c::NO_ENTRY;
        }
        for (Rank rank : bwd_visited_nodes)
        {
            _bwdDistances[rank] = c::NO_ENTRY;
        }
        return best_dist;
    }

    int getDistanceBuckets(Rank source, Rank target, int thread_id)
    {
        vector<Distance> &fwdDistances = _fwd_distances_per_thread[thread_id];
        vector<Distance> &bwdDistances = _bwd_distances_per_thread[thread_id];
        if (source == target)
        {
            return 0;
        }
        vector<vector<Rank>> fwd_buckets(_numHops);
        vector<vector<Rank>> bwd_buckets(_numHops);
        vector<Rank> fwd_visited_nodes;
        vector<Rank> bwd_visited_nodes;
        fwd_visited_nodes.reserve(3000);
        bwd_visited_nodes.reserve(3000);
        Group source_group = _getGroupOfNode(source);
        fwdDistances[source] = 0;
        fwd_visited_nodes.push_back(source);
        for (LabelElement &hub : _fwdLabels[source])
        {
            Distance new_distance = hub.second;
            if (fwdDistances[hub.first] == c::NO_ENTRY)
            {
                fwdDistances[hub.first] = new_distance;
                fwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (source_group < next_group)
                {
                    fwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < fwdDistances[hub.first])
            {
                fwdDistances[hub.first] = new_distance;
            }
        }
        for (Group group = source_group + 1; group < _numHops; group++)
        {
            for (Rank rank : fwd_buckets[group])
            {
                Distance distance = fwdDistances[rank];
                for (LabelElement &hub : _fwdLabels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (fwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        fwdDistances[hub.first] = new_distance;
                        fwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            fwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < fwdDistances[hub.first])
                    {
                        fwdDistances[hub.first] = new_distance;
                    }
                }
            }
        }
        Group target_group = _getGroupOfNode(target);
        bwdDistances[target] = 0;
        bwd_visited_nodes.push_back(target);
        for (LabelElement &hub : _bwdLabels[target])
        {
            Distance new_distance = hub.second;
            if (bwdDistances[hub.first] == c::NO_ENTRY)
            {
                bwdDistances[hub.first] = new_distance;
                bwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (target_group < next_group)
                {
                    bwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < bwdDistances[hub.first])
            {
                bwdDistances[hub.first] = new_distance;
            }
        }
        for (Group group = target_group + 1; group < _numHops; group++)
        {
            for (Rank rank : bwd_buckets[group])
            {
                Distance distance = bwdDistances[rank];
                for (LabelElement &hub : _bwdLabels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (bwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        bwdDistances[hub.first] = new_distance;
                        bwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            bwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < bwdDistances[hub.first])
                    {
                        bwdDistances[hub.first] = new_distance;
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        for (Rank rank : fwd_visited_nodes)
        {
            if (bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = bwdDistances[rank] + fwdDistances[rank];
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                }
            }
            fwdDistances[rank] = c::NO_ENTRY;
        }
        for (Rank rank : bwd_visited_nodes)
        {
            bwdDistances[rank] = c::NO_ENTRY;
        }
        return best_dist;
    }

    Distance getDistanceHubLabels(Rank source, Rank target)
    {
        Distance dist = c::NO_ENTRY;
        int i1 = 0;
        int i2 = 0;
        const vector<LabelElement> &l1 = _fwdLabels[source];
        const vector<LabelElement> &l2 = _bwdLabels[target];
        while (i1 < l1.size() && i2 < l2.size())
        {
            if (l1[i1].first == l2[i2].first)
            {
                Distance new_dist = l1[i1].second + l2[i2].second;
                if (dist == c::NO_ENTRY || new_dist < dist)
                {
                    dist = new_dist;
                }
                i1++;
                i2++;
            }
            else if (l1[i1].first < l2[i2].first)
            {
                i1++;
            }
            else
            {
                i2++;
            }
        }
        return dist;
    }
    vector<LabelElement> getCHPathWithCH(Rank source, Rank target)
    {
        vector<LabelElement> ch_path;
        if (source == target)
        {
            return ch_path;
        }
        vector<vector<Rank>> fwd_buckets(_numHops);
        vector<vector<Rank>> bwd_buckets(_numHops);
        vector<Rank> fwd_visited_nodes;
        vector<Rank> bwd_visited_nodes;
        fwd_visited_nodes.reserve(3000);
        bwd_visited_nodes.reserve(3000);
        Group source_group = _getGroupOfNode(source);
        _fwdDistances[source] = 0;
        fwd_visited_nodes.push_back(source);
        for (LabelElement &hub : _fwdLabels[source])
        {
            Distance new_distance = hub.second;
            if (_fwdDistances[hub.first] == c::NO_ENTRY)
            {
                _fwdDistances[hub.first] = new_distance;
                _added_by_fwd[hub.first] = source;
                fwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (source_group < next_group)
                {
                    fwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < _fwdDistances[hub.first])
            {
                _fwdDistances[hub.first] = new_distance;
                _added_by_fwd[hub.first] = source;
            }
        }
        for (Group group = source_group + 1; group < _numHops; group++)
        {
            for (Rank rank : fwd_buckets[group])
            {
                Distance distance = _fwdDistances[rank];
                for (LabelElement &hub : _fwdLabels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (_fwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        _fwdDistances[hub.first] = new_distance;
                        _added_by_fwd[hub.first] = rank;
                        fwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            fwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < _fwdDistances[hub.first])
                    {
                        _fwdDistances[hub.first] = new_distance;
                        _added_by_fwd[hub.first] = rank;
                    }
                }
            }
        }
        Group target_group = _getGroupOfNode(target);
        _bwdDistances[target] = 0;
        bwd_visited_nodes.push_back(target);
        for (LabelElement &hub : _bwdLabels[target])
        {
            Distance new_distance = hub.second;
            if (_bwdDistances[hub.first] == c::NO_ENTRY)
            {
                _bwdDistances[hub.first] = new_distance;
                _added_by_bwd[hub.first] = target;
                bwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (target_group < next_group)
                {
                    bwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < _bwdDistances[hub.first])
            {
                _bwdDistances[hub.first] = new_distance;
                _added_by_bwd[hub.first] = target;
            }
        }
        for (Group group = target_group + 1; group < _numHops; group++)
        {
            for (Rank rank : bwd_buckets[group])
            {
                Distance distance = _bwdDistances[rank];
                for (LabelElement &hub : _bwdLabels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (_bwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        _bwdDistances[hub.first] = new_distance;
                        _added_by_bwd[hub.first] = rank;
                        bwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            bwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < _bwdDistances[hub.first])
                    {
                        _bwdDistances[hub.first] = new_distance;
                        _added_by_bwd[hub.first] = rank;
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        Rank mid_node = c::NO_ENTRY;
        Distance dist_fwd = c::NO_ENTRY;
        Distance dist_bwd = c::NO_ENTRY;
        for (int i = 0; i < fwd_visited_nodes.size(); i++)
        {
            Rank rank = fwd_visited_nodes[i];
            Distance fwd_distance = _fwdDistances[rank];
            if (_bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = _bwdDistances[rank] + fwd_distance;
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                    mid_node = rank;
                    dist_fwd = fwd_distance;
                    dist_bwd = _bwdDistances[rank];
                }
            }
        }
        if (best_dist != c::NO_ENTRY)
        {
            vector<LabelElement> upward_path;
            upward_path.reserve(100);
            {
                Distance path_dist = 0;
                Rank cur_node = mid_node;
                Distance cur_dist_source = dist_fwd;
                while (_added_by_fwd[cur_node] != c::NO_ENTRY)
                {
                    Rank cur_target = cur_node;
                    cur_node = _added_by_fwd[cur_node];
                    Distance cur_dist = _fwdDistances[cur_node];
                    int edge_id_from_cur_node = c::NO_ENTRY;
                    Distance dist_to_target = _fwdDistances[cur_target] - cur_dist;
                    const vector<LabelElement> &cur_label = _fwdLabels[cur_node];
                    for (int id = cur_label.size() - 1; id >= 1; id--)
                    {
                        const LabelElement &cur_edge = cur_label[id];
                        if (cur_edge.first <= cur_target && cur_edge.second <= dist_to_target && _fwdDistances[cur_edge.first] == cur_dist + cur_edge.second)
                        {
                            if (edge_id_from_cur_node == c::NO_ENTRY && cur_edge.first == cur_target)
                            {
                                edge_id_from_cur_node = id;
                            }
                            else
                            {
                                Distance dist_to_target2 = dist_to_target - cur_edge.second;
                                const vector<LabelElement> &fwd_label = _ch.getNeighbors(cur_edge.first, true);
                                for (const LabelElement &elem : fwd_label)
                                {
                                    if (elem.first == cur_target && elem.second == dist_to_target2)
                                    {
                                        upward_path.push_back(elem);
                                        cur_target = cur_edge.first;
                                        edge_id_from_cur_node = id;
                                        path_dist += elem.second;
                                        dist_to_target -= elem.second;
                                    }
                                }
                            }
                        }
                    }
                    upward_path.push_back(_fwdLabels[cur_node][edge_id_from_cur_node]);
                    path_dist += _fwdLabels[cur_node][edge_id_from_cur_node].second;
                }
            }
            vector<LabelElement> downward_path;
            downward_path.reserve(100);
            {
                Distance path_dist = 0;
                Rank cur_node = mid_node;
                Distance cur_dist_source = dist_bwd;
                while (_added_by_bwd[cur_node] != c::NO_ENTRY)
                {
                    Rank cur_target = cur_node;
                    cur_node = _added_by_bwd[cur_node];
                    Distance cur_dist = _bwdDistances[cur_node];
                    int edge_id_from_cur_node = c::NO_ENTRY;
                    Distance dist_to_target = _bwdDistances[cur_target] - cur_dist;
                    const vector<LabelElement> &cur_label = _bwdLabels[cur_node];
                    for (int id = cur_label.size() - 1; id >= 1; id--)
                    {
                        const LabelElement &cur_edge = cur_label[id];
                        if (cur_edge.first <= cur_target && cur_edge.second <= dist_to_target && _bwdDistances[cur_edge.first] == cur_dist + cur_edge.second)
                        {
                            if (edge_id_from_cur_node == c::NO_ENTRY && cur_edge.first == cur_target)
                            {
                                edge_id_from_cur_node = id;
                            }
                            else
                            {
                                Distance dist_to_target2 = dist_to_target - cur_edge.second;
                                const vector<LabelElement> &bwd_label = _ch.getNeighbors(cur_edge.first, false);
                                for (const LabelElement &elem : bwd_label)
                                {
                                    if (elem.first == cur_target && elem.second == dist_to_target2)
                                    {
                                        downward_path.push_back(elem);
                                        cur_target = cur_edge.first;
                                        edge_id_from_cur_node = id;
                                        path_dist += elem.second;
                                        dist_to_target -= elem.second;
                                    }
                                }
                            }
                        }
                    }
                    downward_path.push_back(_bwdLabels[cur_node][edge_id_from_cur_node]);
                    path_dist += _bwdLabels[cur_node][edge_id_from_cur_node].second;
                }
            }
            ch_path.reserve(upward_path.size() + downward_path.size());
            ch_path.resize(upward_path.size());
            reverse_copy(upward_path.begin(), upward_path.end(), ch_path.begin());
            ch_path.insert(ch_path.end(), downward_path.begin(), downward_path.end());
        }

        for (Rank rank : fwd_visited_nodes)
        {
            _added_by_fwd[rank] = c::NO_ENTRY;
            _fwdDistances[rank] = c::NO_ENTRY;
        }
        for (Rank rank : bwd_visited_nodes)
        {
            _added_by_bwd[rank] = c::NO_ENTRY;
            _bwdDistances[rank] = c::NO_ENTRY;
        }
        return ch_path;
    }

    vector<LabelElement> getCHPath(Rank source, Rank target)
    {
        vector<LabelElement> ch_path;
        if (source == target)
        {
            return ch_path;
        }
        vector<vector<Rank>> fwd_buckets(_numHops);
        vector<vector<Rank>> bwd_buckets(_numHops);
        vector<Rank> fwd_visited_nodes;
        vector<Rank> bwd_visited_nodes;
        fwd_visited_nodes.reserve(3000);
        bwd_visited_nodes.reserve(3000);
        Group source_group = _getGroupOfNode(source);
        _fwdDistances[source] = 0;
        fwd_visited_nodes.push_back(source);
        for (int id = 0; id < _fwdLabels[source].size(); id++)
        {
            LabelElement &hub = _fwdLabels[source][id];
            Distance new_distance = hub.second;
            if (_fwdDistances[hub.first] == c::NO_ENTRY)
            {
                _fwdDistances[hub.first] = new_distance;
                _added_by_fwd_2[hub.first] = make_pair(source, id);
                fwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (source_group < next_group)
                {
                    fwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < _fwdDistances[hub.first])
            {
                _fwdDistances[hub.first] = new_distance;
                _added_by_fwd_2[hub.first] = make_pair(source, id);
            }
        }
        for (Group group = source_group + 1; group < _numHops; group++)
        {
            for (Rank rank : fwd_buckets[group])
            {
                Distance distance = _fwdDistances[rank];
                for (int id = 0; id < _fwdLabels[rank].size(); id++)
                {
                    LabelElement &hub = _fwdLabels[rank][id];
                    Distance new_distance = hub.second + distance;
                    if (_fwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        _fwdDistances[hub.first] = new_distance;
                        _added_by_fwd_2[hub.first] = make_pair(rank, id);
                        fwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            fwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < _fwdDistances[hub.first])
                    {
                        _fwdDistances[hub.first] = new_distance;
                        _added_by_fwd_2[hub.first] = make_pair(rank, id);
                    }
                }
            }
        }
        Group target_group = _getGroupOfNode(target);
        _bwdDistances[target] = 0;
        bwd_visited_nodes.push_back(target);
        for (int id = 0; id < _bwdLabels[target].size(); id++)
        {
            LabelElement &hub = _bwdLabels[target][id];
            Distance new_distance = hub.second;
            if (_bwdDistances[hub.first] == c::NO_ENTRY)
            {
                _bwdDistances[hub.first] = new_distance;
                _added_by_bwd_2[hub.first] = make_pair(target, id);
                bwd_visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (target_group < next_group)
                {
                    bwd_buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < _bwdDistances[hub.first])
            {
                _bwdDistances[hub.first] = new_distance;
                _added_by_bwd_2[hub.first] = make_pair(target, id);
                ;
            }
        }
        for (Group group = target_group + 1; group < _numHops; group++)
        {
            for (Rank rank : bwd_buckets[group])
            {
                Distance distance = _bwdDistances[rank];
                for (int id = 0; id < _bwdLabels[rank].size(); id++)
                {
                    LabelElement &hub = _bwdLabels[rank][id];
                    Distance new_distance = hub.second + distance;
                    if (_bwdDistances[hub.first] == c::NO_ENTRY)
                    {
                        _bwdDistances[hub.first] = new_distance;
                        _added_by_bwd_2[hub.first] = make_pair(rank, id);
                        bwd_visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            bwd_buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < _bwdDistances[hub.first])
                    {
                        _bwdDistances[hub.first] = new_distance;
                        _added_by_bwd_2[hub.first] = make_pair(rank, id);
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        Rank mid_node = c::NO_ENTRY;
        Distance dist_fwd = c::NO_ENTRY;
        Distance dist_bwd = c::NO_ENTRY;
        for (int i = 0; i < fwd_visited_nodes.size(); i++)
        {
            Rank rank = fwd_visited_nodes[i];
            Distance fwd_distance = _fwdDistances[rank];
            if (_bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = _bwdDistances[rank] + fwd_distance;
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                    mid_node = rank;
                    dist_fwd = fwd_distance;
                    dist_bwd = _bwdDistances[rank];
                }
            }
        }
        if (best_dist != c::NO_ENTRY)
        {
            vector<LabelElement> upward_path;
            upward_path.reserve(100);
            {
                Distance path_dist = 0;
                Rank cur_node = mid_node;
                while (_added_by_fwd_2[cur_node].first != c::NO_ENTRY)
                {
                    const pair<Rank, int> &last_ch_edge = _fwd_last_ch_edges[_added_by_fwd_2[cur_node].first][_added_by_fwd_2[cur_node].second];
                    if (last_ch_edge.first == c::NO_ENTRY)
                    {
                        upward_path.push_back(_fwdLabels[_added_by_fwd_2[cur_node].first][_added_by_fwd_2[cur_node].second]);
                        cur_node = _added_by_fwd_2[cur_node].first;
                    }
                    else
                    {
                        upward_path.push_back(_fwdLabels[last_ch_edge.first][last_ch_edge.second]);
                        cur_node = last_ch_edge.first;
                    }
                }
            }
            vector<LabelElement> downward_path;
            downward_path.reserve(100);
            {
                Distance path_dist = 0;
                Rank cur_node = mid_node;
                while (_added_by_bwd_2[cur_node].first != c::NO_ENTRY)
                {
                    const pair<Rank, int> &last_ch_edge = _bwd_last_ch_edges[_added_by_bwd_2[cur_node].first][_added_by_bwd_2[cur_node].second];
                    if (last_ch_edge.first == c::NO_ENTRY)
                    {
                        upward_path.push_back(_bwdLabels[_added_by_bwd_2[cur_node].first][_added_by_bwd_2[cur_node].second]);
                        cur_node = _added_by_bwd_2[cur_node].first;
                    }
                    else
                    {
                        upward_path.push_back(_bwdLabels[last_ch_edge.first][last_ch_edge.second]);
                        cur_node = last_ch_edge.first;
                    }
                }
            }
            ch_path.reserve(upward_path.size() + downward_path.size());
            ch_path.resize(upward_path.size());
            reverse_copy(upward_path.begin(), upward_path.end(), ch_path.begin());
            ch_path.insert(ch_path.end(), downward_path.begin(), downward_path.end());
        }

        for (Rank rank : fwd_visited_nodes)
        {
            _added_by_fwd_2[rank] = make_pair(c::NO_ENTRY, c::NO_ENTRY);
            _fwdDistances[rank] = c::NO_ENTRY;
        }
        for (Rank rank : bwd_visited_nodes)
        {
            _added_by_bwd_2[rank] = make_pair(c::NO_ENTRY, c::NO_ENTRY);
            _bwdDistances[rank] = c::NO_ENTRY;
        }
        return ch_path;
    }

    void switchAddedBy()
    {
        if (_added_by_fwd_2.size() == 0)
        {
            _added_by_fwd.clear();
            _added_by_fwd_2.resize(_numNodes, make_pair(c::NO_ENTRY, c::NO_ENTRY));
            _added_by_bwd.clear();
            _added_by_bwd_2.resize(_numNodes, make_pair(c::NO_ENTRY, c::NO_ENTRY));
        }
        else
        {
            _added_by_fwd_2.clear();
            _added_by_fwd.resize(_numNodes, c::NO_ENTRY);
            _added_by_bwd_2.clear();
            _added_by_bwd.resize(_numNodes, c::NO_ENTRY);
        }
    }

    void writeToBinary(string file)
    {
        cout << "Start binary write" << endl;
        ofstream myHLfile(file, ios_base::out | ios_base::binary);
        boost::archive::binary_oarchive oa(myHLfile);
        oa &(*this);
    }

    double speedTestWithQueryFile(string query_file, bool print = true, int query_type = c::QUERY_BUCKETS)
    {
        Timer timer;
        vector<Rank> sources;
        vector<Rank> targets;
        sources.reserve(1000000);
        targets.reserve(1000000);
        long sumDistances = 0;
        ifstream reader(query_file);
        {
            Rank source, target;
            while (reader >> source >> target)
            {
                sources.push_back(source);
                targets.push_back(target);
            }
        }
        int numRuns = sources.size();
        if (query_type == c::QUERY_BUCKETS)
        {
            if (print)
            {
                cout << "testing getDistanceBuckets" << endl;
            }
            timer.start();
            for (int i = 0; i < numRuns; i++)
            {
                sumDistances += getDistanceBuckets(sources[i], targets[i]);
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_HL)
        {
            if (print)
            {
                cout << "testing getDistanceHubLabels" << endl;
            }
            timer.start();
            for (int i = 0; i < numRuns; i++)
            {
                sumDistances += getDistanceHubLabels(sources[i], targets[i]);
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_CHHL)
        {
            if (print)
            {
                cout << "testing getDistanceCHHL" << endl;
            }
            timer.start();
            for (int i = 0; i < numRuns; i++)
            {
                sumDistances += getDistanceCHHL(sources[i], targets[i]);
            }
            timer.stop();
        }
        if (print)
        {
            cout << "Average query time: " << timer.secs() / numRuns << endl;
        }
        cout << "sum of distances: " << sumDistances << endl;
        return timer.secs() / numRuns;
    }

    double speedTestPathWithQueryFile(string query_file, bool print = true, int query_type = c::NO_ENTRY, int nr_of_threads = c::NO_ENTRY)
    {
        Timer timer;
        vector<Rank> sources;
        vector<Rank> targets;
        sources.reserve(1000000);
        targets.reserve(1000000);
        ifstream reader(query_file);
        {
            Rank source, target;
            while (reader >> source >> target)
            {
                sources.push_back(source);
                targets.push_back(target);
            }
        }
        int numRuns = sources.size();
        long sum_path_lengths = 0;
        if (query_type == c::QUERY_CHPATHWITHCH)
        {
            if (print)
            {
                cout << "testing getCHPath with CH" << endl;
            }
            timer.start();
            for (int i = 0; i < numRuns; i++)
            {
                vector<LabelElement> path = getCHPathWithCH(sources[i], targets[i]);
                sum_path_lengths += path.size();
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_CHPATH)
        {
            cout << "Clear ch auxiliary data..." << endl;
            _ch.clearAuxiliaryData();
            _fwdDistances.clear();
            _fwdDistances.shrink_to_fit();
            _bwdDistances.clear();
            _bwdDistances.shrink_to_fit();
            _added_by_fwd.clear();
            _added_by_fwd.shrink_to_fit();
            _added_by_bwd.clear();
            _added_by_bwd.shrink_to_fit();
            timer.start();
            setLastCHEdges(nr_of_threads);
            timer.stop();
            cout << "Setting ch edges with " << (nr_of_threads == c::NO_ENTRY ? 1 : nr_of_threads) << " threads took " << timer.secs() << " seconds" << endl;
            timer.reset();
            _ch.clear();
            _fwdDistances.resize(_numNodes, c::NO_ENTRY);
            _bwdDistances.resize(_numNodes, c::NO_ENTRY);
            switchAddedBy();
            if (print)
            {
                cout << "testing getCHPath" << endl;
            }
            timer.start();
            for (int i = 0; i < numRuns; i++)
            {
                vector<LabelElement> path = getCHPath(sources[i], targets[i]);
                sum_path_lengths += path.size();
            }
            timer.stop();
            switchAddedBy();
        }

        if (print)
        {
            cout << "Average query time: " << timer.secs() / numRuns << endl;
        }
        cout << "average path length: " << sum_path_lengths / numRuns << endl;
        return timer.secs() / numRuns;
    }
    double getAverageLabelSize(bool print = false)
    {
        vector<long> labels_per_group = getSumLabelSizePerGroup();
        long sum_labels = 0;
        for (long l : labels_per_group)
        {
            sum_labels += l;
        }
        if (print)
        {
            cout << "average label size: " << sum_labels / _numNodes << endl;
        }
        return (double)(sum_labels) / _numNodes;
    }

    double getAverageFwdLabelSize(bool print = false)
    {
        vector<long> labels_per_group = getSumFwdLabelSizePerGroup();
        long sum_labels = 0;
        for (long l : labels_per_group)
        {
            sum_labels += l;
        }
        if (print)
        {
            cout << "average forward label size: " << sum_labels / _numNodes << endl;
        }
        return (double)(sum_labels) / _numNodes;
    }

    vector<long> getSumLabelSizePerGroup()
    {
        vector<long> sums_per_group(_numHops, 0);
        Group group = 0;
        for (long i = 0; i < _numNodes; i++)
        {
            if (i == _groupBounds.at(group))
            {
                group++;
            }

            sums_per_group.at(group) += _fwdLabels.at(i).size();
            sums_per_group.at(group) += _bwdLabels.at(i).size();
        }
        return sums_per_group;
    }

    vector<long> getSumFwdLabelSizePerGroup()
    {
        vector<long> sums_per_group(_numHops, 0);
        Group group = 0;
        for (long i = 0; i < _numNodes; i++)
        {
            if (i == _groupBounds.at(group))
            {
                group++;
            }

            sums_per_group.at(group) += _fwdLabels.at(i).size();
        }
        return sums_per_group;
    }

    void initializeWithCH(string ch_file)
    {
        _ch.readFromBinary(ch_file);
        _numNodes = _ch.nofNodes();
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        _added_by_fwd.assign(_numNodes, c::NO_ENTRY);
        _added_by_bwd.assign(_numNodes, c::NO_ENTRY);
    }

    vector<vector<LabelElement>> getFwdLabels()
    {
        return _fwdLabels;
    }

    vector<vector<LabelElement>> getBwdLabels()
    {
        return _bwdLabels;
    }

    //NOTE: _groupBounds, _rankToNodeID, _nodeIDToRank and _ch must be available
    vector<long> getQueryComplexityWithCH(Rank source, int thread_id = c::NO_ENTRY)
    {
        long num_visited_edges = 0;
        long num_visited_nodes = 1;
        int num_visited_labels = 0;
        long label_size = 0;
        vector<Distance> &fwdDistances = thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread.at(thread_id);
        vector<bool> is_visited(_numNodes, false);
        vector<vector<Rank>> fwd_buckets(_numHops);
        Group source_group = _getGroupOfNode(source);
        fwd_buckets.at(source_group).push_back(source);
        for (Group group = source_group; group < _numHops; group++)
        {
            for (int i = 0; i < fwd_buckets.at(group).size(); i++)
            {
                Rank source_rank = fwd_buckets.at(group).at(i);
                fwdDistances[source_rank] = 0;
                priority_queue<pair<Distance, Rank>, vector<pair<Distance, Rank>>, std::greater<pair<Distance, Rank>>> pq;
                vector<Rank> fwd_visited_nodes;
                fwd_visited_nodes.push_back(source_rank);
                pq.push(make_pair(0, source_rank));
                num_visited_labels++;
                while (!pq.empty())
                {
                    pair<Distance, Rank> element = pq.top();
                    pq.pop();
                    Rank rank = element.second;
                    Distance distance = element.first;
                    if (fwdDistances[rank] != distance)
                    {
                        continue;
                    }
                    const vector<LabelElement> &edges = _ch.getNeighbors(rank, true);
                    for (const LabelElement &edge : edges)
                    {
                        Rank next_rank = edge.first;
                        if (rank < next_rank)
                        {
                            Distance new_disstance = distance + edge.second;
                            if (fwdDistances.at(next_rank) == c::NO_ENTRY || fwdDistances.at(next_rank) > new_disstance)
                            {
                                bool should_be_stalled = false;
                                const vector<LabelElement> &bwd_edges = _ch.getNeighbors(next_rank, false);
                                for (const LabelElement &bwd_edge : bwd_edges)
                                {
                                    if (fwdDistances[bwd_edge.first] != c::NO_ENTRY)
                                    {
                                        if (fwdDistances[bwd_edge.first] + bwd_edge.second < new_disstance)
                                        {
                                            should_be_stalled = true;
                                            break;
                                        }
                                    }
                                }
                                if (!should_be_stalled && _ch.getDistance(source_rank, next_rank, thread_id) == new_disstance)
                                {
                                    Group next_group = _getGroupOfNode(next_rank);
                                    if (fwdDistances[next_rank] == c::NO_ENTRY)
                                    {
                                        fwd_visited_nodes.push_back(next_rank);
                                        if (!is_visited.at(next_rank))
                                        {
                                            num_visited_nodes++;
                                            is_visited.at(next_rank) = true;
                                            if (next_group > group)
                                            {
                                                fwd_buckets[next_group].push_back(next_rank);
                                            }
                                        }
                                    }
                                    fwdDistances[next_rank] = new_disstance;
                                    if (group == next_group)
                                    {
                                        pq.push(make_pair(new_disstance, next_rank));
                                    }
                                }
                            }
                        }
                    }
                }
                num_visited_edges += fwd_visited_nodes.size();
                if (group == source_group)
                {
                    label_size = fwd_visited_nodes.size();
                }
                for (Rank r : fwd_visited_nodes)
                {
                    fwdDistances.at(r) = c::NO_ENTRY;
                }
            }
        }
        vector<long> result;
        result.push_back(num_visited_edges);
        result.push_back(num_visited_nodes);
        result.push_back(num_visited_labels);
        result.push_back(label_size);
        return result;
    }

    //NOTE: _groupBounds, _rankToNodeID, _nodeIDToRank and _ch must be available
    long getLabelSizeWithCH(Rank source, int thread_id = c::NO_ENTRY)
    {
        vector<Distance> &fwdDistances = thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread.at(thread_id);
        Group source_group = _getGroupOfNode(source);

        Rank source_rank = source;
        fwdDistances[source_rank] = 0;
        priority_queue<pair<Distance, Rank>, vector<pair<Distance, Rank>>, std::greater<pair<Distance, Rank>>> pq;
        vector<Rank> fwd_visited_nodes;
        fwd_visited_nodes.push_back(source_rank);
        pq.push(make_pair(0, source_rank));
        while (!pq.empty())
        {
            pair<Distance, Rank> element = pq.top();
            pq.pop();
            Rank rank = element.second;
            Distance distance = element.first;
            if (fwdDistances[rank] != distance)
            {
                continue;
            }
            const vector<LabelElement> &edges = _ch.getNeighbors(rank, true);
            for (const LabelElement &edge : edges)
            {
                Rank next_rank = edge.first;
                if (rank < next_rank)
                {
                    Distance new_disstance = distance + edge.second;
                    if (fwdDistances.at(next_rank) == c::NO_ENTRY || fwdDistances.at(next_rank) > new_disstance)
                    {
                        bool should_be_stalled = false;
                        const vector<LabelElement> &bwd_edges = _ch.getNeighbors(next_rank, false);
                        for (const LabelElement &bwd_edge : bwd_edges)
                        {
                            if (fwdDistances[bwd_edge.first] != c::NO_ENTRY)
                            {
                                if (fwdDistances[bwd_edge.first] + bwd_edge.second < new_disstance)
                                {
                                    should_be_stalled = true;
                                    break;
                                }
                            }
                        }
                        if (!should_be_stalled && _ch.getDistance(source_rank, next_rank, thread_id) == new_disstance)
                        {
                            Group next_group = _getGroupOfNode(next_rank);
                            if (fwdDistances[next_rank] == c::NO_ENTRY)
                            {
                                fwd_visited_nodes.push_back(next_rank);
                            }
                            fwdDistances[next_rank] = new_disstance;
                            if (source_group == next_group)
                            {
                                pq.push(make_pair(new_disstance, next_rank));
                            }
                        }
                    }
                }
            }
        }

        long label_size = fwd_visited_nodes.size();
        for (Rank r : fwd_visited_nodes)
        {
            fwdDistances.at(r) = c::NO_ENTRY;
        }
        return label_size;
    }

    void setGroupBounds(vector<Group> level_bounds)
    {
        _groupBounds.clear();
        for (int i = 0; i < level_bounds.size(); i++)
        {
            _groupBounds.push_back(_ch.getFirstNodeOfLevel(level_bounds[i]));
        }
        _groupBounds.push_back(_numNodes);
        _numHops = _groupBounds.size();
    }

    void setLastCHEdges(int nr_of_threads = c::NO_ENTRY)
    {
        cout << "Start finding last ch edges" << endl;
        if (nr_of_threads == c::NO_ENTRY)
        {
            nr_of_threads = 1;
        }
        _fwd_last_ch_edges.resize(_numNodes);
        _bwd_last_ch_edges.resize(_numNodes);
        Level max_level = _ch.getMaxLevel();
        Rank stop_rank = _numNodes;
        for (Level level = max_level; level >= 0; level--)
        {
            Rank start_rank = _ch.getFirstNodeOfLevel(level);
#pragma omp parallel for num_threads(nr_of_threads)
            for (Rank rank = start_rank; rank < stop_rank; rank++)
            {
                vector<pair<Rank, int>> fwd_last_ch_edges(_fwdLabels[rank].size(), make_pair(c::NO_ENTRY, c::NO_ENTRY));
                {
                    const vector<LabelElement> &ch_fwd_edges = _ch.getNeighbors(rank, true);
                    vector<int> edge_indices(ch_fwd_edges.size(), 1);
                    vector<Distance> neighbor_distances(ch_fwd_edges.size(), c::NO_ENTRY);
                    for (int id = 0; id < _fwdLabels[rank].size(); id++)
                    {
                        const LabelElement &edge = _fwdLabels[rank][id];
                        for (int i = 0; i < ch_fwd_edges.size(); i++)
                        {
                            if (edge.first == ch_fwd_edges.at(i).first && edge.second == ch_fwd_edges.at(i).second)
                            {
                                neighbor_distances.at(i) = edge.second;
                            }
                            if (edge_indices.at(i) < _fwdLabels[ch_fwd_edges[i].first].size())
                            {
                                while (edge_indices.at(i) < _fwdLabels[ch_fwd_edges[i].first].size() && edge.first > _fwdLabels[ch_fwd_edges[i].first].at(edge_indices.at(i)).first)
                                {
                                    edge_indices.at(i)++;
                                }
                                if (edge_indices.at(i) < _fwdLabels[ch_fwd_edges[i].first].size())
                                {
                                    const LabelElement &n_edge = _fwdLabels[ch_fwd_edges[i].first].at(edge_indices.at(i));
                                    if (neighbor_distances.at(i) != c::NO_ENTRY && n_edge.first == edge.first && n_edge.second + neighbor_distances.at(i) == edge.second)
                                    {
                                        if (_fwd_last_ch_edges[ch_fwd_edges[i].first].at(edge_indices.at(i)).first == c::NO_ENTRY)
                                        {
                                            fwd_last_ch_edges.at(id).first = ch_fwd_edges[i].first;
                                            fwd_last_ch_edges.at(id).second = edge_indices.at(i);
                                        }
                                        else
                                        {
                                            fwd_last_ch_edges.at(id) = _fwd_last_ch_edges[ch_fwd_edges[i].first].at(edge_indices.at(i));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                _fwd_last_ch_edges[rank] = fwd_last_ch_edges;
                vector<pair<Rank, int>> bwd_last_ch_edges(_bwdLabels[rank].size(), make_pair(c::NO_ENTRY, c::NO_ENTRY));
                {
                    const vector<LabelElement> &ch_bwd_edges = _ch.getNeighbors(rank, false);
                    vector<int> edge_indices(ch_bwd_edges.size(), 1);
                    vector<Distance> neighbor_distances(ch_bwd_edges.size(), c::NO_ENTRY);
                    for (int id = 0; id < _bwdLabels[rank].size(); id++)
                    {
                        const LabelElement &edge = _bwdLabels[rank][id];
                        for (int i = 0; i < ch_bwd_edges.size(); i++)
                        {
                            if (edge.first == ch_bwd_edges.at(i).first && edge.second == ch_bwd_edges.at(i).second)
                            {
                                neighbor_distances.at(i) = edge.second;
                            }
                            if (edge_indices.at(i) < _bwdLabels[ch_bwd_edges[i].first].size())
                            {
                                while (edge_indices.at(i) < _bwdLabels[ch_bwd_edges[i].first].size() && edge.first > _bwdLabels[ch_bwd_edges[i].first].at(edge_indices.at(i)).first)
                                {
                                    edge_indices.at(i)++;
                                }
                                if (edge_indices.at(i) < _bwdLabels[ch_bwd_edges[i].first].size())
                                {
                                    const LabelElement &n_edge = _bwdLabels[ch_bwd_edges[i].first].at(edge_indices.at(i));
                                    if (neighbor_distances.at(i) != c::NO_ENTRY && n_edge.first == edge.first && n_edge.second + neighbor_distances.at(i) == edge.second)
                                    {
                                        if (_bwd_last_ch_edges[ch_bwd_edges[i].first].at(edge_indices.at(i)).first == c::NO_ENTRY)
                                        {
                                            bwd_last_ch_edges.at(id).first = ch_bwd_edges[i].first;
                                            bwd_last_ch_edges.at(id).second = edge_indices.at(i);
                                        }
                                        else
                                        {
                                            bwd_last_ch_edges.at(id) = _bwd_last_ch_edges[ch_bwd_edges[i].first].at(edge_indices.at(i));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                _bwd_last_ch_edges[rank] = bwd_last_ch_edges;
            }
            stop_rank = start_rank;
        }
    }
private:
    long _numNodes;
    vector<vector<LabelElement>> _fwdLabels;
    vector<vector<LabelElement>> _bwdLabels;
    vector<Distance> _fwdDistances;
    vector<Distance> _bwdDistances;
    vector<Rank> _added_by_fwd;
    vector<Rank> _added_by_bwd;
    vector<pair<Rank, int>> _added_by_fwd_2;
    vector<pair<Rank, int>> _added_by_bwd_2;
    vector<Rank> _groupBounds;
    vector<vector<Distance>> _fwd_distances_per_thread;
    vector<vector<Distance>> _bwd_distances_per_thread;
    int _numHops;
    CH _ch;
    vector<vector<pair<Rank, int>>> _fwd_last_ch_edges;
    vector<vector<pair<Rank, int>>> _bwd_last_ch_edges;

    Group _getGroupOfNode(Rank rank)
    {
        Group group = 0;
        while (rank >= _groupBounds.at(group))
        {
            group++;
        }
        return group;
    }

    void _constructLabels(string ch_binary, vector<Level> group_bounds_in_levels, int nr_of_threads = c::NO_ENTRY, int prune_labels = c::NO_ENTRY)
    {
        cout << "Start constructing labels with " << (nr_of_threads == c::NO_ENTRY ? 1 : nr_of_threads) << " thread(s)." << endl;
        cout << "Levels:";
        for (int i = 0; i < group_bounds_in_levels.size(); i++)
        {
            cout << " " << group_bounds_in_levels[i];
        }
        cout << endl
             << "Pruning: ";
        if (prune_labels == c::NO_ENTRY)
        {
            cout << "full" << endl;
        }
        else if (prune_labels == 0)
        {
            cout << "none" << endl;
            prune_labels = _numHops;
        }
        else
        {
            prune_labels = _numHops - prune_labels;
            cout << "until level " << group_bounds_in_levels.at(prune_labels - 1) << endl;
        }
        Level max_level;
        vector<Rank> first_node_of_level;
        {
            CH ch;
            ch.readFromBinary(ch_binary);
            _numNodes = ch.nofNodes();
            _fwdLabels.resize(_numNodes);
            _bwdLabels.resize(_numNodes);
            for (Rank i = 0; i < _numNodes; i++)
            {
                _fwdLabels[i] = ch.getNeighbors(i, true);
                _bwdLabels[i] = ch.getNeighbors(i, false);
            }
            max_level = ch.getMaxLevel();
            for (int j = 0; j <= max_level; j++)
            {
                first_node_of_level.push_back(ch.getFirstNodeOfLevel(j));
            }
            if (nr_of_threads == c::NO_ENTRY)
            {
                _fwdDistances.assign(_numNodes, c::NO_ENTRY);
                _bwdDistances.assign(_numNodes, c::NO_ENTRY);
            }
        }
        if (nr_of_threads != c::NO_ENTRY)
        {
            _fwdDistances.clear();
            _bwdDistances.clear();
            _fwd_distances_per_thread.resize(nr_of_threads);
            _bwd_distances_per_thread.resize(nr_of_threads);
            for (int i = 0; i < nr_of_threads; i++)
            {
                _fwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
                _bwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
            }
        }
        _groupBounds.clear();
        for (int i = 0; i < group_bounds_in_levels.size(); i++)
        {
            _groupBounds.push_back(first_node_of_level.at(group_bounds_in_levels[i]));
        }
        _groupBounds.push_back(_numNodes);

        Rank i = 0;
        int round = -1;
        int save_prune_labels = prune_labels;
        int group_counter = group_bounds_in_levels.size();
        while (i < _numNodes)
        {
            Rank first_node = i;
            round++;
            if (round == 0 || (group_counter >= 0 && group_bounds_in_levels.at(group_counter) - 1 == max_level - round))
            {
                prune_labels = max_level + 1;
                group_counter--;
            }
            else
            {
                prune_labels = save_prune_labels;
            }
            i = _numNodes - first_node_of_level.at(max_level - round);
            long num_new_shortcuts = 0;
            if (nr_of_threads == c::NO_ENTRY)
            {
                for (Rank x = first_node; x < i; x++)
                {
                    Rank rank = _numNodes - 1 - x;
                    num_new_shortcuts += _constructLabelWithoutCH(rank, c::NO_ENTRY, prune_labels);
                }
            }
            else
            {
#pragma omp parallel for num_threads(nr_of_threads)
                for (Rank x = first_node; x < i; x++)
                {
                    int tid = omp_get_thread_num();
                    Rank rank = _numNodes - 1 - x;
                    num_new_shortcuts += _constructLabelWithoutCH(rank, tid, prune_labels);
                }
            }
        }
    }

    int _constructLabel(Rank rank, int thread_id = c::NO_ENTRY, int prune_labels = c::NO_ENTRY)
    {
        Group group = _getGroupOfNode(rank);
        const vector<LabelElement> &fwdNeighbors = _ch.getNeighbors(rank, true);
        const vector<LabelElement> &bwdNeighbors = _ch.getNeighbors(rank, false);
        int num_edges = fwdNeighbors.size() + bwdNeighbors.size();
        vector<LabelElement> fwdLabel;
        vector<LabelElement> bwdLabel;
        LabelElement defaultLabel;
        defaultLabel.first = rank;
        defaultLabel.second = 0;
        fwdLabel.push_back(defaultLabel);
        bwdLabel.push_back(defaultLabel);
        for (int j = 0; j < fwdNeighbors.size(); j++)
        {
            const LabelElement &edge = fwdNeighbors.at(j);
            Rank neighborRank = edge.first;
            if (neighborRank > rank)
            {
                Group neighbor_group = _getGroupOfNode(neighborRank);
                if (neighbor_group == group)
                {
                    const vector<LabelElement> &l2 = _fwdLabels.at(neighborRank);
                    fwdLabel = _mergeTwoLabels(fwdLabel, l2, edge.second);
                }
                else
                {
                    vector<LabelElement> l2;
                    LabelElement element;
                    element.first = neighborRank;
                    element.second = 0;
                    l2.push_back(element);
                    fwdLabel = _mergeTwoLabels(fwdLabel, l2, edge.second);
                }
            }
        }
        for (int j = 0; j < bwdNeighbors.size(); j++)
        {
            const LabelElement &edge = bwdNeighbors.at(j);
            Rank neighborRank = edge.first;
            if (neighborRank > rank)
            {
                if (neighborRank < _groupBounds.at(group))
                {
                    const vector<LabelElement> &l2 = _bwdLabels.at(neighborRank);
                    bwdLabel = _mergeTwoLabels(bwdLabel, l2, edge.second);
                }
                else
                {
                    vector<LabelElement> l2;
                    LabelElement element;
                    element.first = neighborRank;
                    element.second = 0;
                    l2.push_back(element);
                    bwdLabel = _mergeTwoLabels(bwdLabel, l2, edge.second);
                }
            }
        }
        _fwdLabels.at(rank) = fwdLabel;
        _bwdLabels.at(rank) = bwdLabel;
        if (prune_labels == c::NO_ENTRY || group >= prune_labels)
        {
            _fwdLabels.at(rank) = _pruneLabel(_fwdLabels.at(rank), rank, true, thread_id);
            _bwdLabels.at(rank) = _pruneLabel(_bwdLabels.at(rank), rank, false, thread_id);
        }
        return _fwdLabels.at(rank).size() + _bwdLabels.at(rank).size() - num_edges - 2;
    }

    int _constructLabelWithoutCH(Rank rank, int thread_id = c::NO_ENTRY, int prune_labels = c::NO_ENTRY)
    {
        Group group = _getGroupOfNode(rank);
        int num_edges = _fwdLabels[rank].size() + _bwdLabels[rank].size();
        vector<LabelElement> fwdLabel;
        vector<LabelElement> bwdLabel;
        LabelElement defaultLabel;
        defaultLabel.first = rank;
        defaultLabel.second = 0;
        fwdLabel.push_back(defaultLabel);
        bwdLabel.push_back(defaultLabel);
        for (int j = 0; j < _fwdLabels[rank].size(); j++)
        {
            const LabelElement &edge = _fwdLabels[rank].at(j);
            Rank neighborRank = edge.first;
            if (neighborRank > rank)
            {
                Group neighbor_group = _getGroupOfNode(neighborRank);
                if (neighbor_group == group)
                {
                    const vector<LabelElement> &l2 = _fwdLabels.at(neighborRank);
                    fwdLabel = _mergeTwoLabels(fwdLabel, l2, edge.second);
                }
                else
                {
                    vector<LabelElement> l2;
                    LabelElement element;
                    element.first = neighborRank;
                    element.second = 0;
                    l2.push_back(element);
                    fwdLabel = _mergeTwoLabels(fwdLabel, l2, edge.second);
                }
            }
        }
        for (int j = 0; j < _bwdLabels[rank].size(); j++)
        {
            const LabelElement &edge = _bwdLabels[rank].at(j);
            Rank neighborRank = edge.first;
            if (neighborRank > rank)
            {
                if (neighborRank < _groupBounds.at(group))
                {
                    const vector<LabelElement> &l2 = _bwdLabels.at(neighborRank);
                    bwdLabel = _mergeTwoLabels(bwdLabel, l2, edge.second);
                }
                else
                {
                    vector<LabelElement> l2;
                    LabelElement element;
                    element.first = neighborRank;
                    element.second = 0;
                    l2.push_back(element);
                    bwdLabel = _mergeTwoLabels(bwdLabel, l2, edge.second);
                }
            }
        }
        _fwdLabels.at(rank) = fwdLabel;
        _bwdLabels.at(rank) = bwdLabel;
        if (prune_labels == c::NO_ENTRY || group >= prune_labels)
        {
            _fwdLabels.at(rank) = _pruneLabel(_fwdLabels.at(rank), rank, true, thread_id);
            _bwdLabels.at(rank) = _pruneLabel(_bwdLabels.at(rank), rank, false, thread_id);
        }
        return _fwdLabels.at(rank).size() + _bwdLabels.at(rank).size() - num_edges - 2;
    }

    void _constructLabelHL(Rank rank, int thread_id = c::NO_ENTRY)
    {
        const vector<LabelElement> &fwdNeighbors = _ch.getNeighbors(rank, true);
        const vector<LabelElement> &bwdNeighbors = _ch.getNeighbors(rank, false);
        vector<LabelElement> fwdLabel;
        vector<LabelElement> bwdLabel;
        LabelElement defaultLabel;
        defaultLabel.first = rank;
        defaultLabel.second = 0;
        fwdLabel.push_back(defaultLabel);
        bwdLabel.push_back(defaultLabel);
        for (int j = 0; j < fwdNeighbors.size(); j++)
        {
            const LabelElement &edge = fwdNeighbors.at(j);
            Rank neighborRank = edge.first;
            if (neighborRank > rank)
            {
                const vector<LabelElement> &l2 = _fwdLabels.at(neighborRank);
                fwdLabel = _mergeTwoLabels(fwdLabel, l2, edge.second);
            }
        }
        for (int j = 0; j < bwdNeighbors.size(); j++)
        {
            const LabelElement &edge = bwdNeighbors.at(j);
            Rank neighborRank = edge.first;
            if (neighborRank > rank)
            {
                const vector<LabelElement> &l2 = _bwdLabels.at(neighborRank);
                bwdLabel = _mergeTwoLabels(bwdLabel, l2, edge.second);
            }
        }
        _fwdLabels.at(rank) = fwdLabel;
        _bwdLabels.at(rank) = bwdLabel;
        _fwdLabels.at(rank) = _pruneLabelHL(_fwdLabels.at(rank), rank, true, thread_id);
        _bwdLabels.at(rank) = _pruneLabelHL(_bwdLabels.at(rank), rank, false, thread_id);
    }

    vector<LabelElement>
    _pruneLabel(const vector<LabelElement> &label, Rank rank, bool forward, int thread_id = c::NO_ENTRY)
    {
        vector<LabelElement> prunedLabel;
        prunedLabel.reserve(label.size());
        vector<Rank> visited_nodes = _setDistances(rank, forward, thread_id);
        for (int i = 0; i < label.size(); i++)
        {
            Rank neighborRank = label.at(i).first;
            Distance dist = _getDistanceOneWay(neighborRank, forward, thread_id);
            assert(dist <= label.at(i).second);
            if (dist == label.at(i).second)
            {
                prunedLabel.push_back(label.at(i));
            }
        }
        _resetDistances(visited_nodes, forward, thread_id);

        return prunedLabel;
    }

    vector<LabelElement>
    _pruneLabelHL(const vector<LabelElement> &label, Rank rank, bool forward, int thread_id = c::NO_ENTRY)
    {
        vector<LabelElement> prunedLabel;
        prunedLabel.reserve(label.size());
        for (int i = 0; i < label.size(); i++)
        {
            Rank neighborRank = label.at(i).first;
            if (forward)
            {
                if (getDistanceHubLabels(rank, neighborRank) == label.at(i).second)
                {
                    prunedLabel.push_back(label.at(i));
                }
            }
            else
            {
                if (getDistanceHubLabels(neighborRank, rank) == label.at(i).second)
                {
                    prunedLabel.push_back(label.at(i));
                }
            }
        }

        return prunedLabel;
    }

    vector<LabelElement>
    _mergeTwoLabels(const vector<LabelElement> &l1, const vector<LabelElement> &l2, Distance distance2)
    {
        vector<LabelElement> label;
        label.reserve(l1.size() + l2.size());
        int i1 = 0;
        int i2 = 0;
        while (i2 < l2.size())
        {
            while (i1 < l1.size() && l1[i1].first < l2[i2].first)
            {
                label.push_back(l1[i1++]);
            }
            if (i1 < l1.size() && l1[i1].first == l2[i2].first)
            {
                if (l1[i1].second <= l2[i2].second + distance2)
                {
                    label.push_back(l1[i1]);
                }
                else
                {
                    LabelElement element = l2[i2];
                    element.second += distance2;
                    label.push_back(element);
                }
                i1++;
            }
            else
            {
                LabelElement element = l2[i2];
                element.second += distance2;
                label.push_back(element);
            }
            i2++;
        }
        while (i1 < l1.size())
        {
            label.push_back(l1[i1++]);
        }

        return label;
    }

    vector<Rank> _setDistances(Rank source, bool forward, int thread_id = c::NO_ENTRY)
    {
        vector<Distance> &distances = forward ? (thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread[thread_id]) : (thread_id == c::NO_ENTRY ? _bwdDistances : _bwd_distances_per_thread[thread_id]);
        vector<vector<LabelElement>> &labels = forward ? _fwdLabels : _bwdLabels;
        vector<vector<Rank>> buckets(_numHops);
        vector<Rank> visited_nodes;
        visited_nodes.reserve(3000);
        Group source_group = _getGroupOfNode(source);
        distances[source] = 0;
        visited_nodes.push_back(source);
        for (LabelElement &hub : labels[source])
        {
            Distance new_distance = hub.second;
            if (distances[hub.first] == c::NO_ENTRY)
            {
                distances[hub.first] = new_distance;
                visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (source_group < next_group)
                {
                    buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < distances[hub.first])
            {
                distances[hub.first] = new_distance;
            }
        }
        for (Group group = source_group + 1; group < _numHops; group++)
        {
            for (Rank rank : buckets[group])
            {
                Distance distance = distances[rank];
                for (LabelElement &hub : labels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (distances[hub.first] == c::NO_ENTRY)
                    {
                        distances[hub.first] = new_distance;
                        visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < distances[hub.first])
                    {
                        distances[hub.first] = new_distance;
                    }
                }
            }
        }
        return visited_nodes;
    }

    void _resetDistances(vector<Rank> &visited_nodes, bool forward, int thread_id = c::NO_ENTRY)
    {
        vector<Distance> &distances = forward ? (thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread[thread_id]) : (thread_id == c::NO_ENTRY ? _bwdDistances : _bwd_distances_per_thread[thread_id]);
        for (int i = 0; i < visited_nodes.size(); i++)
        {
            distances[visited_nodes[i]] = c::NO_ENTRY;
        }
    }

    Distance _getDistanceOneWay(Rank target, bool forward, int thread_id = c::NO_ENTRY)
    {
        vector<vector<Rank>> buckets(_numHops);
        vector<Rank> visited_nodes;
        vector<vector<LabelElement>> &labels = !forward ? _fwdLabels : _bwdLabels;
        vector<Distance> &distances = !forward ? (thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread[thread_id]) : (thread_id == c::NO_ENTRY ? _bwdDistances : _bwd_distances_per_thread[thread_id]);
        Group target_group = _getGroupOfNode(target);
        distances[target] = 0;
        visited_nodes.push_back(target);
        for (LabelElement &hub : labels[target])
        {
            Distance new_distance = hub.second;
            if (distances[hub.first] == c::NO_ENTRY)
            {
                distances[hub.first] = new_distance;
                visited_nodes.push_back(hub.first);
                Group next_group = _getGroupOfNode(hub.first);
                if (target_group < next_group)
                {
                    buckets[next_group].push_back(hub.first);
                }
            }
            else if (new_distance < distances[hub.first])
            {
                distances[hub.first] = new_distance;
            }
        }
        for (Group group = target_group + 1; group < _numHops; group++)
        {
            for (Rank rank : buckets[group])
            {
                Distance distance = distances[rank];
                for (LabelElement &hub : labels[rank])
                {
                    Distance new_distance = hub.second + distance;
                    if (distances[hub.first] == c::NO_ENTRY)
                    {
                        distances[hub.first] = new_distance;
                        visited_nodes.push_back(hub.first);
                        Group next_group = _getGroupOfNode(hub.first);
                        if (group < next_group)
                        {
                            buckets[next_group].push_back(hub.first);
                        }
                    }
                    else if (new_distance < distances[hub.first])
                    {
                        distances[hub.first] = new_distance;
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        vector<Distance> &other_distances = forward ? (thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread[thread_id]) : (thread_id == c::NO_ENTRY ? _bwdDistances : _bwd_distances_per_thread[thread_id]);
        for (Rank rank : visited_nodes)
        {
            if (other_distances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = distances[rank] + other_distances[rank];
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                }
            }
            distances[rank] = c::NO_ENTRY;
        }
        return best_dist;
    }

    int binaryEdgeSearch(const vector<LabelElement> &edges, Rank target, Distance distance)
    {
        if (edges.size() == 0)
        {
            return c::NO_ENTRY;
        }
        int index = (edges.size() + 1) / 2;
        int upper_bound = edges.size() - 1;
        int lower_bound = 0;
        bool finished = false;
        while (!finished)
        {
            index = (upper_bound + lower_bound) / 2;
            if (upper_bound <= lower_bound)
            {
                finished = true;
            }
            if (edges[index].first < target)
            {
                lower_bound = index + 1;
            }
            else if (edges[index].first > target)
            {
                upper_bound = index - 1;
            }
            else if (edges[index].second == distance)
            {
                return index;
            }
            else
            {
                return c::NO_ENTRY;
            }
        }
        return c::NO_ENTRY;
    }
};

void readFromBinary(string file, KHopLabeling &khop)
{
    cout << "Start binary read" << endl;
    ifstream myHLfile(file, ios_base::in | ios_base::binary);
    boost::archive::binary_iarchive ia(myHLfile);
    ia &khop;
    cout << "Finished reading." << endl;
}
void writeToBinary(string file, KHopLabeling &khop)
{
    cout << "Start binary write" << endl;
    ofstream myHLfile(file, ios_base::out | ios_base::binary);
    boost::archive::binary_oarchive oa(myHLfile);
    oa &khop;
}

int main(int argc, char *argv[])
{
    uint seed = time(nullptr);
    string input;
    string output;
    int mode = 0;
    int num_hops = 2;
    int num_threads = c::NO_ENTRY;
    int prune_labels = c::NO_ENTRY;
    string queries;
    string ch_file;

    vector<int> level_bounds;
    while (1)
    {
        int result = getopt(argc, argv, "i:o:m:l:q:h:b:s:t:p:c:");
        if (result == -1)
            break; /* end of list */
        switch (result)
        {
        case '?': /* unknown parameter */
            break;
        case ':': /* missing argument of a parameter */
            fprintf(stderr, "missing argument.\n");
            break;
        case 'i':
            input = optarg;
            break;
        case 'o':
            output = optarg;
            break;
        case 'm':
            mode = stoi(optarg);
            break;
        case 'h':
        {
            num_hops = stoi(optarg);
            break;
        }
        case 'q':
            queries = optarg;
            break;
        case 'b':
        {
            stringstream ss(optarg);
            for (int i; ss >> i;)
            {
                level_bounds.push_back(i);
                if (ss.peek() == ',')
                    ss.ignore();
            }
            break;
        }
        case 's':
        {
            seed = stoul(optarg);
            break;
        }
        case 't':
        {
            num_threads = stoi(optarg);
            break;
        }
        case 'c':
        {
            ch_file = optarg;
            break;
        }
        case 'p':
            prune_labels = stoi(optarg);
            break;
        default: /* unknown */
            break;
        }
    }
    srand(seed);
    if (mode == 1)
    {
		if (output.empty())
		{
			output = "results-" + to_string(seed) + ".log";
		}
        KHopLabeling labeling;
        labeling.constructFromCH(input, level_bounds, num_threads, prune_labels);
        labeling.writeToBinary(output);
    }
    else if (mode == 2)
    {
        vector<int> best_bounds;
        KHopLabeling labeling;
        best_bounds = labeling.sampleBestBoundsWithEstimation(input, num_hops, num_threads);
        labeling.constructFromCH(input, best_bounds, num_threads);
        std::size_t stop = input.find_last_of("-");
        string output_file = input.substr(0, stop);
        for (int i = 0; i < best_bounds.size(); i++)
        {
            output_file += "-" + to_string(best_bounds[i]);
        }
        output_file += ".bin";
        labeling.writeToBinary(output_file);
    }
    else if (mode == 3)
    {
        KHopLabeling labeling;
        labeling.initializeWithCH(ch_file);
        labeling.readFromBinary(input);
        double qt1 = labeling.speedTestWithQueryFile(queries, true, c::QUERY_BUCKETS);
        double qt2 = labeling.speedTestPathWithQueryFile(queries, true, c::QUERY_CHPATHWITHCH);
        double qt3 = labeling.speedTestPathWithQueryFile(queries, true, c::QUERY_CHPATH, num_threads);
        std::size_t start = input.find_last_of("/");
        std::size_t stop = input.find(".");
        string output_file = output + "speedtest-" + input.substr(start + 1, stop - start - 1) + ".log";
        cout << "Writing output to " << output_file << endl;
        ofstream writer(output_file);
        writer << "input: " << input << "\n";
        writer << "k = " << labeling.getBounds().size() << "\n";
        writer << "query file: " << queries << "\n";
        writer << "time: " << getTime();
        writer << "average query time distance computation in seconds: " << qt1 << "\n";
        writer << "average query time ch path computation without last ch edge in seconds: " << qt2 << "\n";
        writer << "average query time ch path computation with last ch edge in seconds: " << qt3 << "\n";
        writer << "size: " << labeling.getSizeOfLabeling() / 1024 / 1024 << " MB"
               << "\n";
        writer.close();
    }
    else if (mode == 4)
	{
		if (output.empty())
		{
			std::cerr << "No output path set" << std::endl;
			return -1;
		}
		CH ch;
		ch.constructCH(input);
		ch.writeToBinary(output);
	}

    return 0;
}

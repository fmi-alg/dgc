#ifndef FILE_CH_SEEN
#define FILE_CH_SEEN

#include <iostream>
#include <vector>
#include <random>
#include "c.h"
#include "Timer.h"
#include "CHParser.h"
#include <fstream>
#include <omp.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

using namespace std;

typedef pair<Rank, Distance> LabelElement;

class CH
{
public:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &_fwdLabels;
        ar &_bwdLabels;
        ar &_firstNodeOfLevel;
        ar &_maxLevel;
    }
    CH()
    {
    }

    void constructCH(string file)
    {
        cout << "Loading CH..." << endl;
        CHParser ch_fmi;
        ch_fmi.readFromFMIFile(file);
        cout << "Finished." << endl;
        _numNodes = ch_fmi.nofNodes();
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        vector<NodeID> rankToNodeID = ch_fmi.rankToNodeID;
        _numNodes = rankToNodeID.size();
        _fwdLabels.resize(_numNodes);
        _bwdLabels.resize(_numNodes);
        vector<Rank> nodeIDToRank(_numNodes);
        _maxLevel = ch_fmi.nodeIDToLevel.at(rankToNodeID.at(_numNodes - 1));
        _firstNodeOfLevel.resize(_maxLevel + 1);
        Level currentLevel = 0;
        _firstNodeOfLevel.at(currentLevel) = 0;
        for (Rank i = 0; i < _numNodes; i++)
        {
            NodeID nodeID = rankToNodeID.at(i);
            nodeIDToRank.at(nodeID) = i;
            Level level = ch_fmi.nodeIDToLevel.at(nodeID);
            if (level > currentLevel)
            {
                currentLevel++;
                _firstNodeOfLevel.at(currentLevel) = i;
            }
        }
#pragma omp parallel for
        for (Rank x = 0; x < _numNodes; x++)
        {
            Rank rank = _numNodes - 1 - x;
            NodeID nodeID = rankToNodeID.at(rank);
            vector<EdgeType> fwdNeighbors = ch_fmi.getNeighbors(nodeID, true);
            vector<EdgeType> bwdNeighbors = ch_fmi.getNeighbors(nodeID, false);
            vector<LabelElement> fwdLabel;
            vector<LabelElement> bwdLabel;
            for (int j = 0; j < fwdNeighbors.size(); j++)
            {
                const EdgeType &edge = fwdNeighbors.at(j);
                NodeID neighborID = edge.target;
                Rank neighborRank = nodeIDToRank.at(neighborID);
                if (neighborRank > rank)
                {
                    LabelElement new_hub;
                    new_hub.first = neighborRank;
                    new_hub.second = edge.weight;
                    fwdLabel.push_back(new_hub);
                }
            }
            for (int j = 0; j < bwdNeighbors.size(); j++)
            {
                const EdgeType &edge = bwdNeighbors.at(j);
                NodeID neighborID = edge.source;
                Rank neighborRank = nodeIDToRank.at(neighborID);
                if (neighborRank > rank)
                {
                    LabelElement new_hub;
                    new_hub.first = neighborRank;
                    new_hub.second = edge.weight;
                    bwdLabel.push_back(new_hub);
                }
            }
            sort(fwdLabel.begin(), fwdLabel.end());
            sort(bwdLabel.begin(), bwdLabel.end());
            _fwdLabels.at(rank) = fwdLabel;
            _bwdLabels.at(rank) = bwdLabel;
        }
    }

    void writeOutRankToNodeID(string input, string output)
    {
        cout << "Loading CH..." << endl;
        CHParser ch_fmi;
        ch_fmi.readFromFMIFile(input);
        cout << "Finished." << endl;
        ofstream writer(output);
        for (Rank i = 0; i < ch_fmi.nofNodes(); i++)
        {
            writer << i << " " << ch_fmi.rankToNodeID[i] << "\n";
        }
        writer.close();
    }

    vector<LabelElement> &getNeighbors(Rank rank, bool forward)
    {
        if (forward)
        {
            return _fwdLabels.at(rank);
        }
        return _bwdLabels.at(rank);
    }

    void readFromBinary(string file)
    {
        cout << "Start binary read" << endl;
        ifstream myHLfile(file, ios_base::in | ios_base::binary);
        boost::archive::binary_iarchive ia(myHLfile);
        ia &(*this);
        _numNodes = _fwdLabels.size();
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        _added_by_fwd.assign(_numNodes, c::NO_ENTRY);
        _added_by_bwd.assign(_numNodes, c::NO_ENTRY);
        cout << "Finished reading." << endl;
    }

    void writeToBinary(string file)
    {
        cout << "Start binary write" << endl;
        ofstream myHLfile(file, ios_base::out | ios_base::binary);
        boost::archive::binary_oarchive oa(myHLfile);
        oa &(*this);
    }

    void setNumberOfThreads(int num_threads)
    {
        if (_fwd_distances_per_thread.size() != num_threads || _bwd_distances_per_thread.size() != num_threads)
        {
            _fwd_distances_per_thread.resize(num_threads);
            _bwd_distances_per_thread.resize(num_threads);
            for (int i = 0; i < num_threads; i++)
            {
                _fwd_distances_per_thread.at(i).resize(_numNodes, c::NO_ENTRY);
                _bwd_distances_per_thread.at(i).resize(_numNodes, c::NO_ENTRY);
            }
        }
    }

    int getDistance(int source, int target, int thread_id = c::NO_ENTRY)
    {
        vector<Distance> &fwdDistances = thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread.at(thread_id);
        vector<Distance> &bwdDistances = thread_id == c::NO_ENTRY ? _bwdDistances : _bwd_distances_per_thread.at(thread_id);
        if (source == target)
        {
            return 0;
        }
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
                    for (LabelElement &label : _fwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (fwdDistances[label.first] == c::NO_ENTRY)
                        {
                            fwdDistances[label.first] = new_distance;
                            fwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < fwdDistances[label.first])
                        {
                            fwdDistances[label.first] = new_distance;
                            pQ.push(make_pair(new_distance, label.first));
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
                    for (LabelElement &label : _bwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (bwdDistances[label.first] == c::NO_ENTRY)
                        {
                            bwdDistances[label.first] = new_distance;
                            bwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < bwdDistances[label.first])
                        {
                            bwdDistances[label.first] = new_distance;
                            pQ.push(make_pair(new_distance, label.first));
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

    //DEBUG
    int debugGetDistance(int source, int target, string output_file, int thread_id = c::NO_ENTRY)
    {
        vector<Distance> &fwdDistances = thread_id == c::NO_ENTRY ? _fwdDistances : _fwd_distances_per_thread.at(thread_id);
        vector<Distance> &bwdDistances = thread_id == c::NO_ENTRY ? _bwdDistances : _bwd_distances_per_thread.at(thread_id);
        if (source == target)
        {
            return 0;
        }
        priority_queue<pair<Distance, Rank>, vector<pair<Distance, Rank>>, std::greater<pair<Distance, Rank>>> pQ;
        vector<Rank> fwd_visited_nodes;
        vector<Rank> bwd_visited_nodes;
        vector<Rank> fwd_added_by(_numNodes, c::NO_ENTRY);
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
                    for (LabelElement &label : _fwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (fwdDistances[label.first] == c::NO_ENTRY)
                        {
                            fwdDistances[label.first] = new_distance;
                            fwd_added_by[label.first] = rank;
                            fwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < fwdDistances[label.first])
                        {
                            fwdDistances[label.first] = new_distance;
                            fwd_added_by[label.first] = rank;
                            pQ.push(make_pair(new_distance, label.first));
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
                    for (LabelElement &label : _bwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (bwdDistances[label.first] == c::NO_ENTRY)
                        {
                            bwdDistances[label.first] = new_distance;
                            bwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < bwdDistances[label.first])
                        {
                            bwdDistances[label.first] = new_distance;
                            pQ.push(make_pair(new_distance, label.first));
                        }
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        Rank midpoint = c::NO_ENTRY;

        for (Rank rank : fwd_visited_nodes)
        {
            if (bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = bwdDistances[rank] + fwdDistances[rank];
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                    midpoint = rank;
                }
            }
        }
        vector<LabelElement> path;
        while (midpoint != c::NO_ENTRY)
        {
            path.push_back(make_pair(midpoint, fwdDistances[midpoint]));
            midpoint = fwd_added_by[midpoint];
        }
        ofstream writer(output_file);
        writer << path.size() << "\n";
        for (int i = 0; i < path.size(); i++)
        {
            writer << source << " " << path[i].first << " " << path[i].second << "\n";
        }
        writer.close();
        for (Rank rank : fwd_visited_nodes)
        {
            fwdDistances[rank] = c::NO_ENTRY;
        }
        for (Rank rank : bwd_visited_nodes)
        {
            bwdDistances[rank] = c::NO_ENTRY;
        }

        return best_dist;
    }
    //DEBUG END

    long nofNodes()
    {
        return _numNodes;
    }

    Level getMaxLevel()
    {
        return _maxLevel;
    }

    Rank getFirstNodeOfLevel(Level level)
    {
        return _firstNodeOfLevel.at(level);
    }

    double speedTest(int numRuns = 100000)
    {
        Timer timer;
        srand(time(nullptr));

        vector<int> sources(numRuns);
        vector<int> targets(numRuns);
        vector<int> results(numRuns);
        for (int i = 0; i < numRuns; i++)
        {
            int source = rand() % _numNodes;
            int target = rand() % _numNodes;
            sources[i] = source;
            targets[i] = target;
        }
        timer.start();
        for (int i = 0; i < numRuns; i++)
        {
            results[i] = getDistance(sources[i], targets[i]);
        }
        timer.stop();
        long sumDistances = 0;
        for (int i = 0; i < numRuns; i++)
        {
            sumDistances += results[i];
        }
        cout << "sum of distances: " << sumDistances << endl;
        return timer.secs() / numRuns;
    }

    double speedTestWithQueryFile(string query_file)
    {
        Timer timer;

        vector<int> sources;
        vector<int> targets;
        vector<int> results;
        sources.reserve(1000000);
        targets.reserve(1000000);
        results.reserve(1000000);
        Rank source, target;
        ifstream reader(query_file);
        while (reader >> source >> target)
        {
            sources.push_back(source);
            targets.push_back(target);
        }
        int numRuns = sources.size();
        timer.start();
        for (int i = 0; i < numRuns; i++)
        {
            results[i] = getDistance(sources[i], targets[i]);
        }
        timer.stop();
        long sumDistances = 0;
        for (int i = 0; i < numRuns; i++)
        {
            sumDistances += results[i];
        }
        cout << "sum of distances: " << sumDistances << endl;
        return timer.secs() / numRuns;
    }

    bool debugWithFMI(string file, int numRuns = 1000)
    {
        CHParser ch_fmi;
        ch_fmi.readFromFMIFile(file);
        vector<int> sources(numRuns);
        vector<int> targets(numRuns);
        vector<int> results(numRuns);
        vector<int> results_fmi(numRuns);
        for (int i = 0; i < numRuns; i++)
        {
            int source = rand() % _numNodes;
            int target = rand() % _numNodes;
            sources[i] = source;
            targets[i] = target;
        }
        for (int i = 0; i < numRuns; i++)
        {
            results[i] = getDistance(sources[i], targets[i]);
            results_fmi[i] = ch_fmi.getDistance(ch_fmi.rankToNodeID[sources[i]], ch_fmi.rankToNodeID[targets[i]]);
        }
        for (int i = 0; i < numRuns; i++)
        {
            if (results[i] != results_fmi[i])
            {
                cout << "Error at query " << i + 1 << ": " << sources[i] << " " << targets[i] << " " << results[i] << " " << results_fmi[i] << " " << results[i] - results_fmi[i] << endl;
                return false;
            }
            else
            {
                cout << "Correct: " << results[i] << endl;
            }
        }
        return true;
    }

    vector<LabelElement> getCHPath(Rank source, Rank target)
    {
        vector<LabelElement> ch_path;
        if (source == target)
        {
            return ch_path;
        }
        priority_queue<pair<Distance, Rank>, vector<pair<Distance, Rank>>, std::greater<pair<Distance, Rank>>> pQ;
        vector<Rank> fwd_visited_nodes;
        vector<Rank> bwd_visited_nodes;
        pQ.push(make_pair(0, source));
        _fwdDistances[source] = 0;
        fwd_visited_nodes.push_back(source);
        while (!pQ.empty())
        {
            pair<Distance, Rank> next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (LabelElement &label : _bwdLabels[rank])
                {
                    if (_fwdDistances[label.first] != c::NO_ENTRY)
                    {
                        Distance new_dist = _fwdDistances[label.first] + label.second;
                        if (new_dist < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    for (LabelElement &label : _fwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (_fwdDistances[label.first] == c::NO_ENTRY)
                        {
                            _fwdDistances[label.first] = new_distance;
                            _added_by_fwd[label.first] = rank;
                            fwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < _fwdDistances[label.first])
                        {
                            _fwdDistances[label.first] = new_distance;
                            _added_by_fwd[label.first] = rank;
                            pQ.push(make_pair(new_distance, label.first));
                        }
                    }
                }
            }
        }
        pQ.push(make_pair(0, target));
        _bwdDistances[target] = 0;
        bwd_visited_nodes.push_back(target);
        while (!pQ.empty())
        {
            pair<Distance, Rank> next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_bwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (LabelElement &label : _fwdLabels[rank])
                {
                    if (_bwdDistances[label.first] != c::NO_ENTRY)
                    {
                        Distance new_dist = _bwdDistances[label.first] + label.second;
                        if (new_dist < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    for (LabelElement &label : _bwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (_bwdDistances[label.first] == c::NO_ENTRY)
                        {
                            _bwdDistances[label.first] = new_distance;
                            _added_by_bwd[label.first] = rank;
                            bwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < _bwdDistances[label.first])
                        {
                            _bwdDistances[label.first] = new_distance;
                            _added_by_bwd[label.first] = rank;
                            pQ.push(make_pair(new_distance, label.first));
                        }
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        Rank mid_node = c::NO_ENTRY;
        for (Rank rank : fwd_visited_nodes)
        {
            if (_bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_distance = _bwdDistances[rank] + _fwdDistances[rank];
                if (best_dist == c::NO_ENTRY || new_distance < best_dist)
                {
                    best_dist = new_distance;
                    mid_node = rank;
                }
            }
        }
        vector<LabelElement> upward_path;
        upward_path.reserve(100);
        {
            Rank cur_node = mid_node;
            while (_added_by_fwd[cur_node] != c::NO_ENTRY)
            {
                Rank previous_node = cur_node;
                cur_node = _added_by_fwd[cur_node];
                for (const LabelElement &elem : _fwdLabels[cur_node])
                {
                    if (elem.first == previous_node && _fwdDistances[cur_node] + elem.second == _fwdDistances[previous_node])
                    {
                        upward_path.push_back(elem);
                        break;
                    }
                }
            }
        }
        vector<LabelElement> downward_path;
        downward_path.reserve(100);
        {
            Rank cur_node = mid_node;
            while (_added_by_bwd[cur_node] != c::NO_ENTRY)
            {
                Rank previous_node = cur_node;
                cur_node = _added_by_bwd[cur_node];
                for (const LabelElement &elem : _bwdLabels[cur_node])
                {
                    if (elem.first == previous_node && _bwdDistances[cur_node] + elem.second == _bwdDistances[previous_node])
                    {
                        downward_path.push_back(elem);
                        break;
                    }
                }
            }
        }
        for (Rank rank : fwd_visited_nodes)
        {
            _fwdDistances[rank] = c::NO_ENTRY;
            _added_by_fwd[rank] = c::NO_ENTRY;
        }
        for (Rank rank : bwd_visited_nodes)
        {
            _bwdDistances[rank] = c::NO_ENTRY;
            _added_by_bwd[rank] = c::NO_ENTRY;
        }
        ch_path.resize(upward_path.size());
        reverse_copy(upward_path.begin(), upward_path.end(), ch_path.begin());
        ch_path.insert(ch_path.end(), downward_path.begin(), downward_path.end());
        return ch_path;
    }

    vector<LabelElement> getCHPath(Rank source, Rank target, Rank mid_node, int fwd_index, int bwd_index, const vector<Rank> &fwd_visited_nodes, const vector<Rank> &bwd_visited_nodes, const vector<Distance> &fwd_distances, const vector<Distance> &bwd_distances)
    {
        vector<LabelElement> ch_path;
        vector<LabelElement> upward_path;
        upward_path.reserve(100);
        Distance dist_fwd = fwd_distances[fwd_index];
        Distance dist_bwd = bwd_distances[bwd_index];
        {
            Distance cur_dist = 0;
            Rank cur_node = mid_node;
            Distance cur_dist_source = dist_fwd;
            for (int i = fwd_index - 1; i >= 0; i--)
            {
                if (fwd_distances[i] <= cur_dist_source)
                {
                    Rank rank = fwd_visited_nodes[i];
                    for (const LabelElement &edge : _fwdLabels[rank])
                    {
                        if (edge.first >= cur_node)
                        {
                            if (edge.first == cur_node)
                            {
                                if (fwd_distances[i] + edge.second + cur_dist == dist_fwd)
                                {
                                    upward_path.push_back(edge);
                                    cur_dist += edge.second;
                                    cur_node = rank;
                                    cur_dist_source = fwd_distances[i];
                                }
                            }
                            /*else {
                                break;
                            }*/
                        }
                    }
                }
            }
            //DEBUG
            if (cur_node != source)
            {
                cout << "ERROR fwd: midnode " << mid_node << ", path ";
                for (int i = 0; i < upward_path.size(); i++)
                {
                    cout << upward_path[i].first << ", ";
                }
                cout << endl;
            }
            //DEBUG END
        }
        vector<LabelElement> downward_path;
        downward_path.reserve(100);
        {
            Distance cur_dist = 0;
            Rank cur_node = mid_node;
            Distance cur_dist_target = dist_bwd;
            for (int i = bwd_index - 1; i >= 0; i--)
            {
                if (bwd_distances[i] <= cur_dist_target)
                {
                    Rank rank = bwd_visited_nodes[i];
                    for (const LabelElement &edge : _bwdLabels[rank])
                    {
                        if (edge.first >= cur_node)
                        {
                            if (edge.first == cur_node)
                            {
                                if (bwd_distances[i] + edge.second + cur_dist == dist_bwd)
                                {
                                    downward_path.push_back(edge);
                                    cur_dist += edge.second;
                                    cur_dist_target = bwd_distances[i];
                                    cur_node = rank;
                                }
                            }
                            /*else{
                                break;
                            }*/
                        }
                    }
                }
            }
            //DEBUG
            if (cur_node != target)
            {
                cout << "ERROR bwd: midnode " << mid_node << ", path ";
                for (int i = 0; i < upward_path.size(); i++)
                {
                    cout << upward_path[i].first << ", ";
                }
                cout << endl;
            }
            //DEBUG END
        }
        ch_path.resize(upward_path.size());
        reverse_copy(upward_path.begin(), upward_path.end(), ch_path.begin());
        ch_path.insert(ch_path.end(), downward_path.begin(), downward_path.end());
        return ch_path;
    }

    long nofEdges()
    {
        long sum_fwd = 0;
        for (Rank i = 0; i < _fwdLabels.size(); i++)
        {
            sum_fwd += _fwdLabels[i].size();
        }
        long sum_bwd = 0;
        for (Rank i = 0; i < _bwdLabels.size(); i++)
        {
            sum_bwd += _bwdLabels[i].size();
        }
        return sum_fwd + sum_bwd;
    }

    long nofFwdBwdEdges(bool forward)
    {
        long sum_edges = 0;
        if (forward)
        {
            for (Rank i = 0; i < _fwdLabels.size(); i++)
            {
                sum_edges += _fwdLabels[i].size();
            }
        }
        else
        {
            for (Rank i = 0; i < _bwdLabels.size(); i++)
            {
                sum_edges += _bwdLabels[i].size();
            }
        }
        return sum_edges;
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

    vector<int> analyzeQuery(int numRuns = 10000000, int nr_of_threads = c::NO_ENTRY)
    {
        if (nr_of_threads == c::NO_ENTRY)
        {
            nr_of_threads = 1;
        }
        if (_fwd_distances_per_thread.size() < nr_of_threads)
        {
            _fwd_distances_per_thread.resize(nr_of_threads);
            _bwd_distances_per_thread.resize(nr_of_threads);
            for (int i = 0; i < nr_of_threads; i++)
            {
                _fwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
                _bwd_distances_per_thread[i].resize(_numNodes, c::NO_ENTRY);
            }
        }
        vector<vector<long>> sum_results_per_thread;
        for (int i = 0; i < nr_of_threads; i++)
        {
            vector<long> sum_results(2, 0);
            sum_results_per_thread.push_back(sum_results);
        }
#pragma omp parallel for num_threads(nr_of_threads)
        for (int i = 0; i < numRuns; i++)
        {
            int thread_id = omp_get_thread_num();
            Rank source = rand() % _numNodes;
            Rank target = rand() % _numNodes;
            vector<int> result = _analyzeQuery(source, target, thread_id);
            for (int j = 0; j < result.size(); j++)
            {
                sum_results_per_thread[thread_id].at(j) += result[j];
            }
        }
        vector<long> sum_results_all(2, 0);
        for (int j = 0; j < nr_of_threads; j++)
        {
            sum_results_all[0] += sum_results_per_thread[j][0];
            sum_results_all[1] += sum_results_per_thread[j][1];
        }
        vector<int> average_results(2);
        for (int j = 0; j < sum_results_all.size(); j++)
        {
            average_results[j] = sum_results_all[j] / numRuns;
        }
        return average_results;
    }

    void createDebugFile(string input_file, string output_file)
    {
        ifstream reader(input_file);
        string output = "";
        Rank source, target;
        int counter = 0;
        while (reader >> source >> target)
        {
            Distance distance = getDistance(source, target);
            output += "\n" + to_string(source) + " " + to_string(target) + " " + to_string(distance);
            counter++;
        }
        reader.close();
        output = to_string(counter) + output;
        ofstream writer(output_file);
        writer << output;
        writer.close();
    }

    vector<vector<LabelElement>> getFwdEdges()
    {
        return _fwdLabels;
    }

    vector<vector<LabelElement>> getBwdEdges()
    {
        return _bwdLabels;
    }

    void sortEdges(){
        for(vector<LabelElement>& vec : _fwdLabels){
            sort(vec.begin(), vec.end());
        }
        for(vector<LabelElement>& vec : _bwdLabels){
            sort(vec.begin(), vec.end());
        }
    }

    void clearAuxiliaryData(){
        _fwdDistances.clear();
        _bwdDistances.clear();
        _added_by_bwd.clear();
        _added_by_fwd.clear();
        _fwd_distances_per_thread.clear();
        _bwd_distances_per_thread.clear();
    }

    void clear()
    {
        _fwdLabels.clear();
        _fwdLabels.shrink_to_fit();
        _bwdLabels.clear();
        _bwdLabels.shrink_to_fit();
        _fwdDistances.clear();
        _bwdDistances.clear();
        _added_by_fwd.clear();
        _added_by_bwd.clear();
    }

private:
    long _numNodes;
    Level _maxLevel;
    vector<vector<LabelElement>> _fwdLabels;
    vector<vector<LabelElement>> _bwdLabels;
    vector<Distance> _fwdDistances;
    vector<Distance> _bwdDistances;
    vector<Rank> _added_by_fwd;
    vector<Rank> _added_by_bwd;
    vector<Rank> _firstNodeOfLevel;
    vector<vector<Distance>> _fwd_distances_per_thread;
    vector<vector<Distance>> _bwd_distances_per_thread;

    vector<int> _analyzeQuery(int source, int target, int thread_id = c::NO_ENTRY)
    {
        int num_labels = 0;
        int num_edges = 0;
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
                num_labels += 1;
                num_edges += _bwdLabels[rank].size();
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
                    num_labels += 1;
                    num_edges += _fwdLabels[rank].size();
                    for (LabelElement &label : _fwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (fwdDistances[label.first] == c::NO_ENTRY)
                        {
                            fwdDistances[label.first] = new_distance;
                            fwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < fwdDistances[label.first])
                        {
                            fwdDistances[label.first] = new_distance;
                            pQ.push(make_pair(new_distance, label.first));
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
                num_labels += 1;
                num_edges += _fwdLabels[rank].size();
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
                    num_labels += 1;
                    num_edges += _bwdLabels[rank].size();
                    for (LabelElement &label : _bwdLabels[rank])
                    {
                        Distance new_distance = label.second + distance;
                        if (bwdDistances[label.first] == c::NO_ENTRY)
                        {
                            bwdDistances[label.first] = new_distance;
                            bwd_visited_nodes.push_back(label.first);
                            pQ.push(make_pair(new_distance, label.first));
                        }
                        else if (new_distance < bwdDistances[label.first])
                        {
                            bwdDistances[label.first] = new_distance;
                            pQ.push(make_pair(new_distance, label.first));
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

        vector<int> result;
        result.push_back(num_labels);
        result.push_back(num_edges);
        return result;
    }
};

#endif
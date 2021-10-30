#pragma once

#include <fstream>
#include <vector>
#include <utility>
#include <iostream>
#include <limits>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <queue>
#include <tuple>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "Timer.h"

using namespace std;

typedef double CoordType;
typedef int EdgeCost;

const EdgeCost MAX_EDGE_COST = numeric_limits<EdgeCost>::max();

const NodeID NO_NODE_ID = -1;
const EdgeID NO_EDGE_ID = -1;

class NodeType // not really interesting information
{
public:
	CoordType lat, lon;
	int64_t fmiID;
	int64_t osmID;
	double elev;
	//		string carry;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar &lat;
		ar &lon;
		ar &fmiID;
		ar &osmID;
		ar &elev;
	}
};

class EdgeType
{
public:
	EdgeType()
	{
		source = -1;
		target = -1;
		weight = -1;
		index = -1;
		child_1 = -1;
		child_2 = -1;
	}
	EdgeType(NodeID _src, NodeID _trg, EdgeCost _wght, vector<double> _values, EdgeID _index, EdgeID _child_1 = -1, EdgeID _child_2 = -1)
	{
		source = _src;
		target = _trg;
		weight = _wght;
		index = _index;
		child_1 = _child_1;
		child_2 = _child_2;
	}
	EdgeID index;
	NodeID source, target;
	EdgeCost weight;
	EdgeID child_1, child_2;

	bool operator<(const EdgeType &itm) const
	{
		if (source < itm.source)
			return true;
		else if ((source == itm.source) && (target < itm.target))
			return true;
		else if ((source == itm.source) && (target == itm.target) && (weight < itm.weight))
			return true;
		return false;
	}
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar &source;
		ar &target;
		ar &weight;
	}
};

class CHParser
{
public:
	int numValues;
	vector<NodeType> nodeList;

	vector<EdgeType> edgeList;

	vector<EdgeID> edgeListOut;
	vector<EdgeID> edgeOffsetOut;

	vector<EdgeID> edgeListIn;
	vector<EdgeID> edgeOffsetIn;
	vector<Level> nodeIDToLevel;
	vector<NodeID> rankToNodeID;

	Level maxLevel;

	vector<EdgeType> getNeighbors(NodeID nodeID, bool forward)
	{
		vector<EdgeType> neighbors;
		if (forward)
		{
			for (EdgeID i = edgeOffsetOut.at(nodeID); i < edgeOffsetOut.at(nodeID + 1); i++)
			{
				neighbors.push_back(edgeList.at(i));
			}
		}
		else
		{
			for (EdgeID i = edgeOffsetIn.at(nodeID); i < edgeOffsetIn.at(nodeID + 1); i++)
			{
				neighbors.push_back(edgeList.at(edgeListIn.at(i)));
			}
		}

		return neighbors;
	}

	vector<EdgeType> getEdges(){
		return edgeList;
	}

	Distance getDistance(NodeID source, NodeID target)
	{
		struct PQElement
		{
			NodeID nodeIndex;
			Distance distance;
		};
		struct PQComparer
		{
			bool operator()(const PQElement &lhs, const PQElement &rhs)
			{
				return lhs.distance > rhs.distance;
			}
		};
		priority_queue<PQElement, vector<PQElement>, PQComparer> pQFwd;
		priority_queue<PQElement, vector<PQElement>, PQComparer> pQBwd;
		vector<Distance> distsFwd(nofNodes(), c::NO_ENTRY);
		vector<Distance> distsBwd(nofNodes(), c::NO_ENTRY);
		PQElement sourceElement;
		sourceElement.distance = 0;
		sourceElement.nodeIndex = source;
		pQFwd.push(sourceElement);
		while (!pQFwd.empty())
		{
			PQElement elem = pQFwd.top();
			pQFwd.pop();
			if (distsFwd[elem.nodeIndex] == c::NO_ENTRY)
			{
				distsFwd[elem.nodeIndex] = elem.distance;
				vector<EdgeType> neighbors = getNeighbors(elem.nodeIndex, true);
				for (int i = 0; i < neighbors.size(); i++)
				{
					EdgeType &edge = neighbors.at(i);
					if (distsFwd[edge.target] == c::NO_ENTRY && nodeIDToLevel[edge.target] > nodeIDToLevel[elem.nodeIndex])
					{
						PQElement new_elem;
						new_elem.distance = edge.weight + elem.distance;
						new_elem.nodeIndex = edge.target;
						pQFwd.push(new_elem);
					}
				}
			}
		}

		PQElement targetElement;
		targetElement.distance = 0;
		targetElement.nodeIndex = target;
		pQBwd.push(targetElement);
		Distance best_dist = c::NO_ENTRY;
		while (!pQBwd.empty())
		{
			PQElement elem = pQBwd.top();
			pQBwd.pop();
			if (distsBwd[elem.nodeIndex] == c::NO_ENTRY)
			{
				distsBwd[elem.nodeIndex] = elem.distance;
				if (distsFwd[elem.nodeIndex] != c::NO_ENTRY && (best_dist == c::NO_ENTRY || distsFwd[elem.nodeIndex] + elem.distance < best_dist)){
					best_dist = distsFwd[elem.nodeIndex] + elem.distance;
				}
				vector<EdgeType> neighbors = getNeighbors(elem.nodeIndex, false);
				for (int i = 0; i < neighbors.size(); i++)
				{
					EdgeType &edge = neighbors.at(i);
					if (distsBwd[edge.source] == c::NO_ENTRY && nodeIDToLevel[edge.source] > nodeIDToLevel[elem.nodeIndex])
					{
						PQElement new_elem;
						new_elem.distance = edge.weight + elem.distance;
						new_elem.nodeIndex = edge.source;
						pQBwd.push(new_elem);
					}
				}
			}
		}

		return best_dist;
	}

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar &nodeList;
		ar &edgeList;
		ar &edgeListOut;
		ar &edgeOffsetOut;
		ar &edgeListIn;
		ar &edgeOffsetIn;
	}

	long nofNodes() const
	{
		return nodeList.size();
	}
	long nofEdges() const
	{
		return edgeListOut.size();
	}

	void readFromFMIFile(string fname)
	{
		Timer myTimer;
		myTimer.start();

		nodeList.clear();

		edgeList.clear();

		edgeListOut.clear();
		edgeOffsetOut.clear();

		edgeListIn.clear();
		edgeOffsetIn.clear();

		ifstream inputFile(fname.c_str());
		char junkC[1024];
		int junkI;
		string line;
		while(true){
			std::getline(inputFile, line);
			if (line.find('#') == std::string::npos && line.size()){
				break;
			}
		}
	
		long nofNodes, nofEdges;
		nofNodes = stoi(line);
		inputFile >> nofEdges;

		cout << "We will read " << nofNodes << " nodes and " << nofEdges << " edges" << endl;

		cout << "Resident memory req should be: " << (nofNodes * (sizeof(NodeType)) +					 // node stuff
													  nofEdges * (sizeof(EdgeType) + sizeof(EdgeType)) + // 2 EdgeLists
													  nofNodes * 2 * sizeof(NodeID) +						 // 2 offset arrays
													  nofEdges * 2 * sizeof(EdgeID)						 // edge lists
													  ) / (1024 * 1024)
			 << "MB" << endl;
		cout << "NodeType: " << sizeof(NodeType) << endl;
		cout << "EdgeType: " << sizeof(EdgeType) << endl;
		cout << "long	: " << sizeof(long) << endl;
		cout << "int	: " << sizeof(int) << endl;
		NodeType curNode;

		nodeList.resize(nofNodes);
		nodeIDToLevel.resize(nofNodes);
		rankToNodeID.resize(nofNodes);
		maxLevel = 0;

		cout << endl
			 << "NODES: ";
		for (long i = 0; i < nofNodes; i++)
		{
			inputFile >> curNode.fmiID;
			inputFile >> curNode.osmID;
			inputFile >> curNode.lat;
			inputFile >> curNode.lon;
			inputFile >> curNode.elev;
			Level level;
			inputFile >> level;
			nodeIDToLevel.at(i) = level;
			if (level > maxLevel)
			{
				maxLevel = level;
			}

			// cout<<"Read: "<<curNodeExt.fmiID<<"\n";
			if ((i) % max(nofNodes / 10, (long) 1) == 1)
				cout << int((i * 100.0) / nofNodes) << "% " << flush;
			inputFile.getline(junkC, 256);
			//				curNode.carry=string(junkC);
			nodeList[i] = curNode;
			rankToNodeID.at(i) = i;
		}
		//compute ranks
		sort(rankToNodeID.begin(), rankToNodeID.end(), [&](NodeID i, NodeID j) { return nodeIDToLevel[i] < nodeIDToLevel[j]; });

		edgeOffsetOut.resize(nofNodes + 1);
		edgeOffsetOut[0] = 0;

		vector<int> nofIncoming(nofNodes, 0);

		NodeID lastSrcID = -1;
		int countZeroWeight = 0;
		int countShortcuts = 0;
		EdgeType curEdge;
		curEdge.weight = 0;
		cout << endl
			 << "EDGES: ";
		for (long j = 0; j < nofEdges; j++)
		{
			NodeID srcID;
			inputFile >> srcID;
			curEdge.source = srcID;
			inputFile >> curEdge.target;
			inputFile >> curEdge.weight;
			int tmp;
			inputFile >> tmp;
			inputFile >> tmp;
			inputFile >> curEdge.child_1;
			inputFile >> curEdge.child_2;
			if ((j) % max(nofEdges / 10, (long) 1) == 1)
				cout << int((j * 100.0) / nofEdges) << "% " << flush;

#ifdef DEBUG2
			cout << "The rest is: " << junkC << endl;

			cout << "Adding Edge from " << srcID
				 << " to " << curEdge.target << " with cost " << curEdge.weight << endl;
#endif
			assert(srcID >= lastSrcID);
			lastSrcID = srcID;
			curEdge.index = j;
			edgeList.push_back(curEdge);
			edgeListOut.push_back(j);
			edgeOffsetOut[srcID + 1] = j + 1;
			nofIncoming[curEdge.target]++;
		}
		edgeOffsetOut[nodeList.size()] = edgeListOut.size();

		cout << endl
			 << "We augmented " << countZeroWeight << " zero-weight edges" << endl;
		cout << endl
			 << "We have " << countShortcuts << " many Shortcuts and hence " << nofEdges - countShortcuts << " original edges" << endl;
#ifdef DEBUG2
		cout << "NodeList ist jetzt " << nodeList.size() << " groß" << endl;
#endif
		cout << "After first Offset creation" << endl;

		// if node v has no outgoing edges, edgeOffsetOut[v+1] is
		// never properly set
		// we need to set edgeOffsetOut[v+1] to edgeOffsetOut[v'] where
		// v'<v+1 and v'-1 had outgoing edges
		for (long i = 0; i < nofNodes; i++)
			if (edgeOffsetOut[i + 1] == 0)
				edgeOffsetOut[i + 1] = edgeOffsetOut[i];

		// we finished construction of outEdge lists and offsets

		// now offsets and stuff for incoming edge list
		edgeOffsetIn.resize(nofNodes + 1, 0);
		edgeListIn.resize(nofEdges);
		edgeOffsetIn[0] = 0;
		for (long i = 0; i < nofNodes; i++)
			edgeOffsetIn[i + 1] = edgeOffsetIn[i] + nofIncoming[i];
		assert(edgeOffsetIn[nofNodes] == nofEdges);

		// now store incoming edges; need to iterate over nodes to get source
		for (long i = 0; i < nofNodes; i++)
		{
			for (EdgeID j = edgeOffsetOut[i]; j < edgeOffsetOut[i + 1]; j++)
			{
				curEdge = edgeList[edgeListOut[j]];
				assert(curEdge.source == i);
				NodeID src = curEdge.source;
				NodeID trg = curEdge.target;

				assert(nofIncoming[trg] > 0);

				edgeListIn[edgeOffsetIn[trg + 1] - nofIncoming[trg]] = edgeListOut[j];
				nofIncoming[trg] = nofIncoming[trg] - 1;
			}
		}
#ifdef DEBUG2
		for (int i = 0; i < nofNodes + 1; i++)
			cout << "offsetInc[" << i << "]=" << edgeOffsetIn[i] << endl;
#endif

		inputFile.close();
		myTimer.stop();
		cout << "Read graph in " << myTimer.secs() << "s" << endl;
	}
};

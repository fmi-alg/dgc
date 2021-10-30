#ifndef FILE_C_SEEN
#define FILE_C_SEEN

#include <climits>

namespace c
{
    int NO_ENTRY = -1;
    int QUERY_MERGE = 0;
    int QUERY_PQ = 1;
    int QUERY_BUCKETS = 2;
    int QUERY_HL = 3;
    int QUERY_CH = 4;
    int QUERY_CHHL = 5;
    int QUERY_BUCKETSMAPS = 6;
    int QUERY_PQSoD = 7;
    int QUERY_CHPATH = 8;
    int QUERY_CHPATHWITHCH = 9;
}

typedef int Rank;
typedef int NodeID;
typedef long EdgeID;
typedef int Distance;
typedef int Level;
typedef int Group;

#endif
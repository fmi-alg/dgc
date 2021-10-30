# Docker Image for DGC

DGC can also be run using a docker container.
This makes it easy to run the full pipeline starting with an URI to a osm.pbf file from OpenStreetMap and ending with running the benchmarks.
Let's give it a try!

## Building

You can build the image or fetch it from Docker Hub.
Note that the latest image on Docker Hub may not reflect the state of the repository.

```bash
docker-compose build
```

This may take quite some time.
You can speed-up things by using the image from Docker Hub:

```bash
docker-compose pull
```

## Configuration

By default all files are stored below the `data` folder.
You may change this default by editing the `.env` file or copy it to another location.
This is especially helpful if you want to benchmark multiple data sets.

```bash
mkdir data/lichtenstein
cp .env lichtenstein.env
$ vim lichtenstein.env
DGC_BASE="./data/lichtenstein"
DGC_SOURCE="${DGC_BASE}/source"
DGC_DATA="${DGC_BASE}/data"
DGC_RESULTS="${DGC_BASE}/results"
DGC_DB="${DGC_BASE}/db"
```

## Import data

To import data you can either place an osm.pbf file into the source folder an rename it to source.osm.pbf or let the docker container take care of downloading a file:

```bash
$ docker-compose --env-file lichtenstein.env run dgc fetch "http://download.geofabrik.de/europe/liechtenstein-latest.osm.pbf"
--2021-10-30 19:39:41--  http://download.geofabrik.de/europe/liechtenstein-latest.osm.pbf
Resolving download.geofabrik.de (download.geofabrik.de)... 95.216.115.119, 116.202.112.212
Connecting to download.geofabrik.de (download.geofabrik.de)|95.216.115.119|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 2499183 (2.4M) [application/octet-stream]
Saving to: '/source/source.osm.pbf'

/source/source.osm.pbf                          100%[====================================================================================================>]   2.38M  6.48MB/s    in 0.4s    

2021-10-30 19:39:41 (6.48 MB/s) - '/source/source.osm.pbf' saved [2499183/2499183]
```

## Extract a Graph

The following command will extract the largest connected component found in the `source.osm.pbf` file.

```bash
$ docker-compose --env-file lichtenstein.env run dgc graph
Extracting largest components
Found 16 config entries


Calculating min/max node id for direct hash map: 0 seconds for 2M 499K 8 
Min nodeId=26860698
Max nodeId=9207797074
There are not enough nodes in the data set to warrant the usage of a direct mapped cache


Collecting candidate node refs: 0 seconds for 2M 499K 8 

Finding unavailable nodes: 2346933|2499008=93.91% (1|0|1)
Finding unavailable nodes: 1 seconds for 2M 499K 8 =2499008 1/s


Collecting needed node refs: 0 seconds for 2M 499K 8 


Collecting nodes: 0 seconds for 24K 214 
Graph has 24214 nodes and 48563 edges.


Writing out nodes: 0 seconds for 24K 214 


Processing ways: 0 seconds for 2M 499K 8 
Finding connected components for 24214 nodes and 48563 edges
Creating single sets for union find
Uniting all nodes using edges
Setting node representatives
Found 20 connected components
Sorting 24214 nodes according to their connected component
Sorting 48563 edges according to their connected component
Found 1 connected components above your threshold


Writing connected components: 0 seconds for 72K 777
```

## Compute a Contraction Hierarchy

We now have to compute a contraction hierarchy:

```bash
$ docker-compose --env-file lichtenstein.env run dgc ch
Computing contraction hierarchy
Removed 6 duplicate edge(s) and updated edge IDs.
Converting CH into binary format
Loading CH...
We will read 24040 nodes and 86552 edges
Resident memory req should be: 9MB
NodeType: 40
EdgeType: 40
long    : 8
int     : 4

NODES: 0% 10% 20% 30% 40% 50% 60% 70% 80% 90% 
EDGES: 0% 10% 20% 30% 40% 50% 59% 69% 79% 89% 99% 
We augmented 0 zero-weight edges

We have 0 many Shortcuts and hence 86552 original edges
After first Offset creation
Read graph in 0.120927s
Finished.
Start binary write
```

## Compute the DGC

The following command will compute the distance closures.
The command needs the ch levels used to partition the ch graph.
They to be given as a single string separated by a comma.

```bash
$ docker-compose --env-file lichtenstein.env run dgc dgc "5,10,15"
Computing dgc with bounds 5,10,15
Start constructing labels with 4 thread(s).
Levels: 5 10 15
Pruning: full
Start binary read
Finished reading.
Finished construction in 0.412014 seconds
Start binary write
```

## Compute Out-of-Memory Data Structures

In order to run the out-of-memory benchmarks we first have to compute the necessary files.
This is also needed to compute sample queries.

```bash
$ docker-compose --env-file lichtenstein.env run dgc oom
Converting to out-of-memory structures
Operation: export-oa
File: /data/ch.bin
dca: hop
Importing data
Finished importing data
Exporting data
Serializing LabelBuckets at 0
Serializing LabelBuckets at 442394
Finished exporting data
Operation: export-oa
File: /data/dgc.bin
dca: hop
Importing data
Finished importing data
Exporting data
Serializing LabelBuckets at 0
Serializing LabelBuckets at 1039650
Finished exporting data
```

## Run the Benchmarks

We're no ready to run the benchmarks.
We first have to generate the queries (you may also provide one yourself):

```bash
$ docker-compose --env-file lichtenstein.env run dgc queries --query-sizes 3 4 5
Setting query_sizes to [3, 4, 5]
Writing files to directory /results/queries/
Generating 10^3 queries
Generating 10^4 queries
Generating 10^5 queries
```

We can now run the benchmarks.
The following command will run the dgc and ch out-of-memory benchmarks each twice with queries of size $10^3$ and $10^4$.
After each iteration the current result table is printed which we omit here for brevity.

```bash
$ docker-compose --env-file lichtenstein.env run dgc bench --runs 2 --query-sizes 3 4 --selection dgc-oom ch-oom

```
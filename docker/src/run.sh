#!/bin/bash
if [[ $EUID -ne 0 ]]; then
   echo "This script must be run as root" 
   exit 1
fi

source /run.env.sh

sourceDir="/source"
dataDir="/data"
resultsDir="/results"
dgcPgPassword="${PGPASSWORD:-dgc}"

# set -x


function createPostgresConfig() {
  cp /etc/postgresql/current/main/postgresql.custom.conf.tmpl /etc/postgresql/current/main/conf.d/postgresql.custom.conf || exit 1
  if [ -n "$AUTOVACUUM" ]; then
    sudo -u postgres echo "autovacuum = $AUTOVACUUM" >> /etc/postgresql/current/main/conf.d/postgresql.custom.conf || exit 1
  fi
  # cat /etc/postgresql/current/main/conf.d/postgresql.custom.conf
}

function setPostgresPassword() {
    sudo -u postgres psql -c "CREATE USER dgc" > /dev/null 2>&1
    sudo -u postgres psql -c "ALTER USER dgc PASSWORD '${dgcPgPassword}'" > /dev/null 2>&1
    sudo -u postgres psql -c "ALTER USER dgc CREATEDB LOGIN" > /dev/null 2>&1
    sudo -u postgres psql -c "GRANT ALL ON TABLESPACE pg_default TO dgc" > /dev/null 2>&1
}

function die() {
    echo ${1}
    service postgresql stop
    exit 1
}

if [ "$1" = "-h" ]; then
    echo "usage: <clean|fetch|graph|ch|dgc|oom|bench>"
    echo "commands:"
    echo "    clean: Clean data directory and erase data base"
    echo "    fetch <uri>: Fetch osm.pbf from the given uri. Defaults to Andorra from Geofabrik"
    echo "    graph: Extract largest connected component from OpenStreetMap file"
    echo "    ch: Create contraction hierarchy based on graph extracted from OpenStreetMap file"
    echo "    dgc <bounds>: Create DGC based on CH data. Bounds defined by a comma-seperated string e.g. \"3,10,15,20\""
    echo "    oom: Create OOM files based on DGC data"
    echo "    create <bounds>: Do steps ch, dhc, oom"
    echo "    queries --query-sizes <nums>: Generate queries with 10^n_i queries"
    echo "    bench [--query-sizes <nums> --runs <num=1>] [--select dgc-oom dgc-im ch-oom ch-im sqloa sql plsql] [-v]: Run selected benchmarks"
    exit 1
fi

fix_permissions.sh || die "Could not fix permissions"

do_graph="n"
do_ch="n"
do_dgc="n"
do_oom="n"
do_queries="n"
do_bench="n"


if [ "$1" = "clean" ]; then
    echo "Cleaning data base related files"
    sudo -u dgc rm  ${dataDir}/* > /dev/null 2>&1
    sudo -u postgres rm -rf /var/lib/postgresql/* > /dev/null 2>&1
fi

if [ "$1" = "fetch" ]; then
    sudo -u dgc wget -O ${sourceDir}/source.osm.pbf "${2}"
    exit 0
fi
if [ "$1" = "graph" ]; then
    do_graph="y"
fi
if [ "$1" = "ch" ]; then
    do_ch="y"
fi
if [ "$1" = "dgc" ]; then
    do_dgc="y"
fi
if [ "$1" = "oom" ]; then
    do_oom="y"
fi
if [ "$1" = "queries" ]; then
    do_queries="y"
fi
if [ "$1" = "bench" ]; then
    do_bench="y"
fi
if [ "$1" = "create" ]; then
    if [ -z "$2" ]; then
        echo "No bounds given"
        exit 1
    fi
    do_graph="y"
    do_ch="y"
    do_dgc="y"
    do_oom="y"
fi

if [ "$do_graph" = "y" ]; then
    #Get the graph
    echo "Extracting largest components"
    sudo -u dgc graph-creator -g fmimaxspeedtext -t time -hs auto -cc topk 1 -c /etc/graph-creator/configs/car.cfg -o ${dataDir}/ ${sourceDir}/source.osm.pbf || die "Failed to extract graph"
    sudo -u dgc mv ${dataDir}/0.cc ${dataDir}/graph.txt || die "Could not rename graph file"
fi

if [ "$do_ch" = "y" ]; then
    echo "Computing contraction hierarchy"
    sudo -u dgc ch-constructor -i ${dataDir}/graph.txt -f FMI -o ${dataDir}/ch.txt -g FMI_CH -t ${CH_CONSTRUCTOR_NUM_THREADS:-4} || die "Failed to compute contraction hierarchy"
    echo "Converting CH into binary format"
    sudo -u dgc dgc-create -m 4 -i ${dataDir}/ch.txt -o ${dataDir}/ch.bin || die "Failed to convert ch graph to binary format"
fi

if [ "$do_dgc" = "y" ]; then
    if [ -z "$2" ]; then
        echo "No bounds given"
        exit 1
    fi
    echo "Computing dgc with bounds $2"
    sudo -u dgc dgc-create -m 1 -i ${dataDir}/ch.bin -o ${dataDir}/dgc.bin -b $2 -t ${DGC_CREATE_NUM_THREADS:-4} || die "Failed to compute dgc"
fi

if [ "$do_oom" = "y" ]; then
    echo "Converting to out-of-memory structures"
    sudo -u dgc dgc-query export-oa -i ${dataDir}/ch.bin || die "Failed to convert ch structures to out-of-memory structures"
    sudo -u dgc dgc-query export-oa -i ${dataDir}/dgc.bin || die "Failed to convert dgc structures to out-of-memory structures"
fi

query_sizes=()
bench_selection=()
runs=1

# Parse remaining command line flags
# Current $1 points to the command
shift
while (( "$#" )); do
  case "$1" in
    -q|--query-sizes)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        #$1 points to the option
        shift
        while (( "$#" )); do
            if [ ${1:0:1} = "-" ]; then
                break
            else
                query_sizes+=($1)
                shift
            fi
        done
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -s|--select|--selection)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        #$1 points to the option
        shift
        while (( "$#" )); do
            if [ ${1:0:1} = "-" ]; then
                break
            else
                bench_selection+=($1)
                shift
            fi
        done
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -r|--runs)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        runs=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -v*)
        verbose_mode="$1"
        shift
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      shift
      ;;
  esac
done

if [ "$do_queries" = "y" ]; then
    mkdir -p ${resultsDir}/queries || die "Could not create queries directory"
    sudo -u dgc bench.py \
        -e dgc-query \
        -q ${resultsDir}/queries/ --generate --query-sizes ${query_sizes[@]} \
        --dgc-oom ${dataDir}/dgc.bin.sserialize-oa -v || die "Generating queries failed"
fi

if [ "$do_bench" = "y" ]; then
    if [[ "${bench_selection[@]}" =~ sql ]]; then
        # Ensure that database directory is in right state
        rm -rf /var/lib/postgresql/* > /dev/null 2>&1
        chown postgres:postgres -R /var/lib/postgresql
        if [ ! -f /var/lib/postgresql/${POSTGRES_VERSION}/main/PG_VERSION ]; then
            sudo -u postgres /usr/lib/postgresql/${POSTGRES_VERSION}/bin/pg_ctl -D /var/lib/postgresql/${POSTGRES_VERSION}/main/ initdb -o "--locale C.UTF-8" > /dev/null 2>&1
        fi
        # Initialize PostgreSQL
        createPostgresConfig
        service postgresql start
        sudo -u postgres createuser dgc > /dev/null 2>&1
        sudo -u postgres createdb -E UTF8 -O dgc dgc > /dev/null 2>&1
        setPostgresPassword
    fi

    mkdir -p ${resultsDir}/bench || die "Could not create bench directory"
    echo "Running benchmarks..."
    sudo -u dgc bench.py \
        -q ${resultsDir}/queries/ --query-sizes ${query_sizes[@]} \
        -e dgc-query \
        --bench ${bench_selection[@]} --runs $runs \
        --drop-caches-cmd /usr/local/bin/drop_caches.sh \
        --checkpoint-dir ${resultsDir}/bench/ -e dgc-query \
        --checkpoint-size 100 \
        --dbuser dgc --dbpass "${dgcPgPassword}" --dbimport --dbdrop \
        --dgc-oom ${dataDir}/dgc.bin.sserialize-oa \
        --dgc-im ${dataDir}/dgc.bin \
        --ch-oom ${dataDir}/ch.bin.sserialize-oa \
        --ch-im ${dataDir}/ch.bin \
        ${verbose_mode} \
        || die "Benchmarking failed"


    if [[ "${bench_selection[@]}" =~ sql ]]; then
        service postgresql stop
    fi
fi
exit 0

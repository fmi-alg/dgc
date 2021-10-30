#!/usr/bin/env python3
# kate: space-indent on; indent-width 4; mixedindent off; indent-mode python;
import argparse
import psycopg2
import re
import sys
import logging
from functools import reduce

cmdLineParser = argparse.ArgumentParser(description="Creates a sql function to compute shortest path routes")
cmdLineParser.add_argument('-d', help='Drop function', dest='drop', action='store_true')
cmdLineParser.add_argument('-g', help='Number of groups', dest='num_groups', type=int)
cmdLineParser.add_argument('--plsql', help='Use pl/pgsql', dest='plsql', action='store_true')
cmdLineParser.add_argument('-q', help='Query source target', dest='query', nargs=2, type=int)
cmdLineParser.add_argument('--query-file', help='Query file', dest='query_file', nargs=1, type=str)
cmdLineParser.add_argument('-e', help='Explain query', dest='explain', action='store_true')
cmdLineParser.add_argument('-b', help='Bench query. Set to stdout to print to stdout.', dest='bench', type=str)
cmdLineParser.add_argument('--batch', help='Process queries in batches of size', dest='batch', type=int)
cmdLineParser.add_argument('--dbname', help='Upload to data base', dest='dbname', type=str)
cmdLineParser.add_argument('--dbuser', help='Database password', dest='dbuser', type=str)
cmdLineParser.add_argument('--dbpass', help='Database password', dest='dbpass', type=str)
cmdLineParser.add_argument('--dbhost', help='Database password', dest='dbhost', type=str, default="localhost")
cmdLineParser.add_argument('--dbport', help='Database password', dest='dbport', type=int, default=5432)
cmdLineParser.add_argument('--quiet', help='Quiet', dest='quiet', action='store_true')
cmdLineParser.add_argument('-v', help='Verbosity', dest='verbosity', action='count', default=0)

parsedArgs = cmdLineParser.parse_args()
num_groups : int

if parsedArgs.verbosity > 1:
    logging.basicConfig(level=logging.DEBUG, format="%(message)s")
elif parsedArgs.verbosity == 1:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
else:
    logging.basicConfig(level=logging.WARN, format="%(message)s")

def connectdb():
    conn_str = "host={host} port={port} user={user} password={password} dbname={dbname} connect_timeout=10".format(
        host=parsedArgs.dbhost,
        port=parsedArgs.dbport,
        user=parsedArgs.dbuser,
        password=parsedArgs.dbpass,
        dbname=parsedArgs.dbname)
    return psycopg2.connect(conn_str)

if parsedArgs.num_groups is None:
    if parsedArgs.dbuser is None:
        logging.error("You need to specify the number of groups");
        sys.exit(1)
    logging.info("Retrieving number of groups from data base")
    with connectdb() as conn:
        conn.set_session(autocommit=True, isolation_level='READ UNCOMMITTED')
        with conn.cursor() as cur:
            cur.execute("select stats.group_count from stats limit 1;")
            result = cur.fetchall()
            num_groups = result[0][0]
else:
    num_groups = parsedArgs.num_groups

if not parsedArgs.quiet:
    logging.info("Batch size: {}".format(parsedArgs.batch))

drop="""drop FUNCTION shop_dist(integer,integer);"""

start="""
CREATE OR REPLACE FUNCTION shop_dist(IN in_src_node int, IN in_tgt_node int)
--    RETURNS table (
--        nid int,
--        dst int
--    )
    returns int
    AS
    $$
DECLARE
    src_group int := 0;
    tgt_group int := 0;
    return_value int := 0;
begin"""

grp_select = """
    select
        g.grp
    into
        src_group
    from
        grp as g
    where 
    (g.node = in_src_node)
    limit 1;

    select
        g.grp
    into
        tgt_group
    from
        grp as g
    where 
    (g.node = in_tgt_node)
    limit 1;"""

tmp_table_def="""
    CREATE TEMP TABLE fwd
    (
        tgt_id int not null,
        dst int not null
    ) ON COMMIT DROP;
    CREATE TEMP TABLE bwd
    (
        tgt_id int not null,
        dst int not null
    ) ON COMMIT DROP;"""

tmp_idx="""
   create unique index fwd_idx on fwd using btree(tgt_id asc nulls last);
   create unique index bwd_idx on bwd using btree(tgt_id asc nulls last);
"""

initial_hopfwd="""
    insert INTO
        fwd(tgt_id, dst)
    select
        fw.tgt, fw.dst 
    from
        fwdlabels{grp} as fw
    where
        fw.src = in_src_node
    ;"""
initial_hopbwd="""
    insert INTO
        bwd(tgt_id, dst)
    select
        bw.tgt, bw.dst
    from
        bwdlabels{grp} as bw
    where
        bw.src = in_tgt_node
    order by bw.tgt
    ;"""

def initial_hop(grp):
    result=[]
    if parsedArgs.plsql:
        result.append("if src_group = {} then".format(grp))
    result.append(initial_hopfwd.format(grp=grp))
    if parsedArgs.plsql:
        result.append("end if;")
        result.append("if tgt_group = {} then".format(grp))
    result.append(initial_hopbwd.format(grp=grp))
    if parsedArgs.plsql:
        result.append("end if;")
    return "\n".join(result)

hopfwd="""
    insert into
        fwd(tgt_id, dst)
    select
        fw.tgt, min(fwd.dst+fw.dst)
    from
        fwd,
        fwdlabels{grp} as fw
    where
        (fwd.tgt_id = fw.src)
    group by (fw.tgt)
    on conflict (tgt_id)
    do update set
        dst = LEAST(fwd.dst, excluded.dst) 
    ;"""
hopbwd="""
    insert into
        bwd(tgt_id, dst)
    select
        bw.tgt, min(bwd.dst+bw.dst)
    from
        bwd,
        bwdlabels{grp} as bw
    where
        (bwd.tgt_id = bw.src)
    group by (bw.tgt)
    on conflict (tgt_id)
    do update set
        dst = LEAST(bwd.dst, excluded.dst) 
    ;"""

def hop(grp):
    result=[]
    if parsedArgs.plsql:
        result.append("if src_group < {} then".format(grp))
    result.append(hopfwd.format(grp=grp))
    if parsedArgs.plsql:
        result.append("end if;")
        result.append("if tgt_group < {} then".format(grp))
    result.append(hopbwd.format(grp=grp))
    if parsedArgs.plsql:
        result.append("end if;")
    return "\n".join(result)


final="""
    --return query 
        select
            min(fwd.dst+bwd.dst)
        into return_value
        from
            fwd,
            bwd
        where
            fwd.tgt_id = bwd.tgt_id
        limit 1
        ;
    drop table fwd;
    drop table bwd;
    return return_value;
END
$$ LANGUAGE {language};"""

explain="""
LOAD 'auto_explain';
SET auto_explain.log_min_duration = 1; -- exclude very fast trivial queries
SET auto_explain.log_nested_statements = ON; -- statements inside functions
SET auto_explain.log_analyze = ON; -- get actual times, too"""

query="""select * from shop_dist({src}, {tgt});"""


def gen_func_str():
    result = []
    if parsedArgs.drop:
        result.append(drop)

    result.append(start)
    # print(grp_select)
    result.append(tmp_table_def)
    #The following does not work since we have to create the table index AFTER the initial step
    #On 1239-1-3-8-30:
    #If we create the tmp index before: 340ms
    #If we do it afterwards: 13ms

    #We expand the table in every hop
    #This has the following advantage:
    #Suppose the group of src node is 3
    #if we first did the inital hop and the expand all the way from group 0
    #We would do work on groups 0-2 for nothing, since no matching nodes are in the table
    #However if we first insert all from group i then the following holds:
    #If src is in group i then all further expansion steps are needed
    #If the group of src is larger than the current group, then the expansion step is cheap since the distance table is empty
    #If the group of src is smaller than the current group, then the insertion step is cheap since that group does not contain src
    #This way we don't need control statements and everything is plain SQL
    if parsedArgs.plsql:
        result.append(grp_select)
    for i in range(0, num_groups):
        result.append(initial_hop(i))
    result.append(tmp_idx)
    for i in range(1, num_groups):
        result.append(hop(i))
    result.append(final.format(language="plpgsql" if parsedArgs.plsql else "plpgsql"))
    return "\n".join(result)

def gen_query():
    result = []
    result.append(query.format(src=parsedArgs.query[0], tgt=parsedArgs.query[1]))
    return "\n".join(result)

def queries_from_file():
    result = []
    with open(parsedArgs.query_file[0]) as f:
        for line in f:
            line : str
            src, tgt = line.split()
            result.append((int(src), int(tgt)))
    logging.debug("Parsed {} queries from file {}".format(len(result), parsedArgs.query_file[0]))
    return result

func_str=gen_func_str()

#Result table contains tuple of (count, time) pairs
result_table = []

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

bench_result_re = re.compile(r'.*Execution\ Time:\s*(\d+\.\d+)')
def parse_bench_result(text: str):
    matches = bench_result_re.match(text)
    if matches.group(1) is None:
        logging.error("No timing information:\n{}".format(text))
        sys.exit(1)
    return int(float(matches.group(1))*1000)

def issue_query(queries: [(int, int)]):
    query_str = "DROP TABLE IF EXISTS queries_table;\n"
    query_str += "CREATE TEMP TABLE queries_table (src int, tgt int) ON COMMIT DROP;\n"
    query_str += "INSERT INTO queries_table (src, tgt) VALUES {};\n".format(
        ",".join(
            ["({}, {})".format(x[0], x[1]) for x in queries]
        )
    )
    if parsedArgs.explain:
        query_str += "explain analyze\n"
    elif parsedArgs.bench:
        query_str += "EXPLAIN (ANALYZE, COSTS OFF, TIMING OFF)\n"
    # query_str += "select queries_table.src, queries_table.tgt, shop_dist(queries_table.src, queries_table.tgt) from queries_table;\n"
    query_str += "select sum(shop_dist(queries_table.src, queries_table.tgt)) from queries_table;\n"
    logging.debug(query_str)
    cur.execute(query_str)
    return str(cur.fetchall())

if parsedArgs.dbname and len(parsedArgs.dbname):
    with connectdb() as conn:
        conn.set_session(autocommit=True, isolation_level='READ UNCOMMITTED')
        with conn.cursor() as cur:
            logging.debug(func_str)
            logging.debug("Uploading function")
            cur.execute(func_str)
            if parsedArgs.explain:
                cur.execute(explain) 
            if parsedArgs.query:
                query_str = gen_query()
                if parsedArgs.explain:
                    query_str = "explain analyze\n{}".format(query_str)
                elif parsedArgs.bench:
                    query_str = "EXPLAIN (ANALYZE, COSTS OFF, TIMING OFF)\n{}".format(query_str)
                logging.debug(query_str)
                result = cur.execute(query_str)
                logging.info(cur.fetchall())
            if parsedArgs.query_file:
                if parsedArgs.batch:
                    queries = list(chunks(queries_from_file(), parsedArgs.batch))
                    num_queries_computed = 0
                    for next_batch in queries:
                        result = issue_query(next_batch)
                        num_queries_computed += len(next_batch)
                        logging.debug(result)
                        if parsedArgs.bench:
                            result_table.append((num_queries_computed, parse_bench_result(result)))

                    cum_time: int = 0
                    for i in range(0, len(result_table)):
                        result_table[i] = (result_table[i][0], result_table[i][1] + cum_time)
                        cum_time = result_table[i][1]
                    result_text = ["sum of distances: -1"]
                    result_text += ["mean time per query [us]: {}".format(int(result_table[-1][1]/result_table[-1][0]))]
                    result_text += ["count,cumulative time"]
                    result_text += ["{}, {}".format(x[0], x[1]) for x in result_table]
                    if parsedArgs.bench == "stdout":
                        print("\n".join(result_text))
                    else:
                        with open(parsedArgs.bench, mode='wt') as f:
                            f.write("\n".join(result_text))
                else:
                    queries=[query.format(x[0], x[1]) for x in queries_from_file()]
                    for q in queries:
                        cur.execute(q)
                        logging.debug(cur.fetachall())
else:
    print(func_str)
    if parsedArgs.query:
        if parsedArgs.explain:
            result.append("explain analyze")
        print(gen_query())

#!/usr/bin/env python3

import typing
from typeguard import typechecked
import logging
import os
import subprocess
import sys
import argparse
import random
import itertools
import psycopg2
from prettytable import PrettyTable
from functools import reduce

query_sizes = [3, 4, 5, 6]

random.seed(1337)

bench_algos = ["dgc-oom", "chhl-oom", "ch-oom", "dgc-im", "chhl-im", "ch-im", "sql", "plsql", "sqloa"]

cmdLineParser = argparse.ArgumentParser(description="Creates benchmark queries for all used data sets. Use @file_name to import command line arguments from file", fromfile_prefix_chars='@')
cmdLineParser.add_argument('-q', help='Queries directory', dest='dest', type=str, required=True)
cmdLineParser.add_argument('--generate', help='Generate queries based on file given by --dgc', dest='gen', action='store_true')
cmdLineParser.add_argument('--bench', help='Benchmark', dest='bench', nargs='+', type=str, choices=["all"] + bench_algos, default=[])
cmdLineParser.add_argument('--dio', help='Enable directio benchmarks', dest='dio', action='store_true')
cmdLineParser.add_argument('-e', help='Path to dgc-query executable', dest='exe', type=str, required=True)
cmdLineParser.add_argument('--dgc-oom', help='Path to out-of-memory dgc data', dest='dgc_oom', type=str, required=True)
cmdLineParser.add_argument('--chhl-oom', help='Path to out-of-memory chhl data', dest='chhl_oom', type=str)
cmdLineParser.add_argument('--ch-oom', help='Path to out-of-memory ch data', dest='ch_oom', type=str)
cmdLineParser.add_argument('--dgc-im', help='Path to dgc in-memory data archive', dest='dgc_im', type=str)
cmdLineParser.add_argument('--chhl-im', help='Path to chhl in-memory data archive', dest='chhl_im', type=str)
cmdLineParser.add_argument('--ch-im', help='Path to ch in-memory data archive', dest='ch_im', type=str)
cmdLineParser.add_argument('--drop-caches-cmd', help='Command to drop cache to bench cold-cache', dest='dropcache', type=str)
cmdLineParser.add_argument('--query-sizes', help='Size of queries in 10^k', dest='query_sizes', nargs='+', type=int)
cmdLineParser.add_argument('--checkpoint-size', help='Measure time after this many queries', dest='checkpoint', type=int)
cmdLineParser.add_argument('--runs', help='How often a query set should be repeated', dest='runs', type=int, default=1)
cmdLineParser.add_argument('--checkpoint-dir', help='Directory where checkpoint timings are stored', dest='checkpoint_dir', type=str)
cmdLineParser.add_argument('--dbname', help='Database name. Defaults to the basename of the file given by --dgc', dest='dbname', type=str)
cmdLineParser.add_argument('--dbuser', help='Database user', dest='dbuser', type=str, default="dgc")
cmdLineParser.add_argument('--dbpass', help='Database password', dest='dbpass', type=str)
cmdLineParser.add_argument('--dbhost', help='Database host', dest='dbhost', type=str, default="localhost")
cmdLineParser.add_argument('--dbport', help='Database port', dest='dbport', type=int, default=5432)
cmdLineParser.add_argument('--dbspace', help='Database tablespace name', dest='dbspace', type=str, default="pg_default")
cmdLineParser.add_argument('--dbimport', help='Import database during benchmarkig. Drops existing database', dest='dbimport', action='store_true')
cmdLineParser.add_argument('--dbdrop', help='Drop database after benchmarkig', dest='dbdrop', action='store_true')
cmdLineParser.add_argument('-v', help='Verbosity', dest='verbosity', action='count', default=0)
args = cmdLineParser.parse_args()

if args.verbosity >= 1:
	logging.basicConfig(level=logging.DEBUG, format="%(message)s")
else:
	logging.basicConfig(level=logging.INFO, format="%(message)s")

if args.bench and "all" in args.bench:
	args.bench = bench_algos

if args.query_sizes:
	logging.debug("Setting query_sizes to {}".format(args.query_sizes))
	query_sizes = args.query_sizes


if args.exe:
	if args.checkpoint:
		if args.checkpoint_dir is None:
			logging.error("You requested checkpoints but did not specify a checkpoint directory")
			sys.exit(1)
		elif not os.path.isdir(args.checkpoint_dir):
			logging.error("Error: checkpoint-dir is not a directory: {}".format(args.checkpoint_dir))
			sys.exit(1)


infiles = {}

if args.dgc_oom:
	infiles["dgc-oom"] = os.path.abspath(args.dgc_oom)
if args.chhl_oom:
	infiles["chhl-oom"] = os.path.abspath(args.chhl_oom)
if args.ch_oom:
	infiles["ch-oom"] = os.path.abspath(args.ch_oom)
if args.dgc_im:
	infiles["dgc-im"] = os.path.abspath(args.dgc_im)
if args.chhl_im:
	infiles["chhl-im"] = os.path.abspath(args.chhl_im)
if args.ch_im:
	infiles["ch-im"] = os.path.abspath(args.ch_im)

logging.debug("Writing files to directory {}".format(args.dest))

@typechecked
def query_fn(num_queries: int):
	return os.path.join(args.dest, "10E{}.txt".format(num_queries))

if args.gen:
	if not os.path.isdir(args.dest):
		logging.error("Destination is not a directory")
		sys.exit(-1)
	
	num_nodes: int = 0
	
	result = subprocess.run([args.exe, "oam", "-i", infiles["dgc-oom"], "--stats"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	result.check_returncode()
	result = result.stdout.split('\n')
	for line in result:
		line: str
		if line.startswith("nodes:"):
			num_nodes = int(line[len("nodes:"):])
			break

	if num_nodes == 0:
		logging.error("Could not determine number of nodes of the graph")
		sys.exit(1)

	for num_queries in query_sizes:
		logging.info("Generating 10^{} queries".format(num_queries))
		with open(query_fn(num_queries), mode="wt") as f:
			for i in range(0, 10**num_queries):
				src = random.randint(0, num_nodes)
				tgt = random.randint(0, num_nodes)
				f.write("{}\t{}\n".format(src, tgt))
				

class CommandExecutionException(Exception):
	pass


def drop_caches():
	logging.debug("Dropping caches")
	result = subprocess.run(["sudo", args.dropcache], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	if result.returncode != 0:
		raise CommandExecutionException("Failed to drop caches with error {}".format(result.stdout))


class Results:
	
	def __init__(self):
		self.tbl = PrettyTable()
		self.tbl.field_names = ["algo", "#q 10^", "async", "dio", "ad rnd", "distance", "time"]
		self.tbl.align = "r"
		self.tbl.align["async"] = "c"
		self.tbl.align["dio"] = "c"
		self.tbl.align["ad rnd"] = "c"
		self.run = None
		self.results = None

	@typechecked
	def b2x(self, v: bool):
		return "X" if v else ""

	@typechecked
	def begin_run(self, algo: str, num_queries: int, threads: int, dio: bool, adv_rnd: bool):
		if not self.run is None:
			raise Exception("begin_run called with a run already in progress")
		self.run = [algo, num_queries, self.b2x(threads > 1), self.b2x(dio), self.b2x(adv_rnd)]
		self.results = []
		self.checkpoint_fn = "{algo}_{num_queries}_{threads}_{dio}_{advise_random}".format(
			algo=algo, num_queries=num_queries, threads=threads, dio=dio, advise_random=adv_rnd)
	
	@typechecked
	def add_result(self, dist: int, time: int, checkpoints):
		if args.checkpoint and checkpoints is not None:
			fn = os.path.join(args.checkpoint_dir, "{}_{}".format(self.checkpoint_fn, len(self.results)))
			with open(fn, mode='wt') as checkpoint_f:
				checkpoint_f.writelines(["{}\n".format(x) for x in checkpoints])
		self.results.append((dist, time))

	def end_run(self):
		dist : int
		time : int
		(dist, time) = reduce(lambda x, y: (x[0]+y[0], x[1]+y[1]), self.results)
		dist /= len(self.results)
		time /= len(self.results)
		if dist != self.results[0][0]:
			logging.error("Distance differs between runs")
		self.run += [dist, time]
		self.tbl.add_row(self.run)
		self.run = None
		self.results = None
		print(self.tbl)

	def write_summary(self):
		if args.checkpoint_dir:
			with open(os.path.join(args.checkpoint_dir, "summary.txt"), 'wt') as summary_fn:
				summary_fn.write("{}".format(self.tbl))

results = Results()


class Bench:
	@typechecked
	def run_process(self, cmd: typing.List[str]) -> typing.List[str]:
		logging.debug("Executing {}".format(" ".join(cmd)))
		with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) as process:
			cmd_output = []
			try:
				for line in process.stdout:
					line = line.rstrip()
					logging.debug(line)
					cmd_output.append(line)
			except:
				process.kill()
				raise
			process.wait()
			if process.returncode != 0:
				raise CommandExecutionException(
					"{output}\n{stderr}\nCommand execution failed with return value of {retval}.\n{cmdline}".format(
						retval=process.returncode, cmdline=" ".join(cmd), stderr=process.stderr.readlines(), output="\n".join(cmd_output)))
			return cmd_output

	#You have to call results.begin_run before
	@typechecked
	def run(self, cmd: typing.List[str]) -> None:
		for run in range(0, args.runs):
			if args.dropcache:
				drop_caches()
			logging.debug(" ".join(cmd))
			cmd_output = self.run_process(cmd)
			cmd_results = []
			for num, line in enumerate(cmd_output):
				num: int
				line: str
				if line.startswith("sum of distances: "):
					cmd_results += [int(line[len("sum of distances: "):])]
				elif line.startswith("mean time per query [us]: "):
					cmd_results += [int(line[len("mean time per query [us]: "):])]
					cmd_results += cmd_output[num+1:]
					break
			results.add_result(int(cmd_results[0]), int(cmd_results[1]), cmd_results[2:] if len(cmd_results) > 2 else [])
		results.end_run()


class DBBench(Bench):
	
	@staticmethod
	@typechecked
	def dbname_from_filename(filepath: str, oa: bool):
		dbname : str = os.path.basename(filepath)
		dbname = dbname[:dbname.find('.')]
		if oa:
			dbname += "-oa"
		else:
			dbname += "-bkt"
		return dbname

	@typechecked
	def __init__(self, name: str, user: str, pw: str, host: str, port: int, tablespace: str):
		self.host = host
		self.port = port
		self.name = name
		self.user = user
		self.pw = pw
		self.tablespace = tablespace

	def connect(self):
		conn = psycopg2.connect(user=self.user, password=self.pw, host=self.host, port=self.port)
		conn.autocommit = True
		return conn

	def dropdb(self):
		logging.debug("Droping database {}".format(self.name))
		conn = None
		try: #https://stackoverflow.com/a/68112827
			conn = self.connect()
			with conn.cursor() as cur:
				cur.execute("DROP DATABASE IF EXISTS \"{}\"".format(self.name))
		finally:
			if conn:
				conn.close()

	@typechecked
	def createdb(self, tablespace: str = "pg_default"):
		logging.debug("Creating database {}".format(self.name))
		conn = None
		try: #https://stackoverflow.com/a/68112827
			conn = self.connect()
			with conn.cursor() as cur:
				cur.execute("""CREATE DATABASE \"{name}\" WITH 
								OWNER = {owner}
								ENCODING = 'UTF8'
								LC_COLLATE = 'C.UTF-8'
								LC_CTYPE = 'C.UTF-8'
								TABLESPACE = {space}
								CONNECTION LIMIT = -1;""".format(name=self.name, owner=self.user, space=tablespace))
		finally:
			if conn:
				conn.close()

	@typechecked
	def _import_data(self, filepath: str, export_cmd: str):
		self.dropdb()
		self.createdb(self.tablespace)
		cmd = [
			args.exe,
			export_cmd,
			"--dbname",
			self.name,
			"--dbhost",
			self.host,
			"--dbport",
			str(self.port),
			"--dbuser",
			self.user,
			"--dbpass",
			self.pw,
			"-i",
			filepath
		]
		cmd_result = super().run_process(cmd)
	
	@typechecked
	def import_data_bucketed(self, filepath: str):
		self._import_data(filepath, "export-sql")
	
	@typechecked
	def import_data_oa(self, filepath: str):
		self._import_data(filepath, "export-sql-oa")
	
	@typechecked
	def run_oa(self, num_queries: int) -> None:
		#sqloa -i ${DGC_DATA} --dbname ${GRAPH}-oa --dbuser dgc --dbpass testtesttest --query-file ${QUERY_FILE} --checkpoint-size ${BATCH_SIZE} --quiet
		cmd = [
			args.exe,
			"sqloa",
			"--dbname",
			self.name,
			"--dbhost",
			self.host,
			"--dbport",
			str(self.port),
			"--dbuser",
			self.user,
			"--dbpass",
			self.pw,
			"-i",
			args.dgc_oom,
			"--query-file",
			query_fn(num_queries),
			"--checkpoint-size",
			str(args.checkpoint),
			"--quiet"
		]
		if args.checkpoint:
			cmd += ["--checkpoint-size", "{}".format(args.checkpoint)]
		
		results.begin_run("sqloa", num_queries, 1, False, False)
		super().run(cmd)

	@typechecked
	def run_sql(self, num_queries: int, plsql: bool) -> None:
		#./sql_shopa_func.py -g ${GROUP_COUNT} --batch ${BATCH_SIZE} -b ${PLSQL_OUT}/10E${QUERY_COUNT}_$i.txt --dbname ${GRAPH}-${SEPARATORS} --query-file ${QUERY_FILE} --plsql
		
		cmd = [
			os.path.join(os.path.dirname(os.path.realpath(__file__)), "sql_shopa_func.py"),
			"--batch",
			str(args.checkpoint),
			"-b",
			"stdout",
			"--dbname",
			self.name,
			"--dbhost",
			self.host,
			"--dbport",
			str(self.port),
			"--dbuser",
			self.user,
			"--dbpass",
			self.pw,
			"--query-file",
			query_fn(num_queries),
			"--quiet"
		]
		algo = "sql"
		if plsql:
			cmd += ["--plsql"]
			algo = "plsql"
		results.begin_run(algo, num_queries, 1, False, False)
		super().run(cmd)


class OOMBench(Bench):

	@typechecked
	def construct_cmdline(self, exe: os.path, infile: os.path, algo: str, queries: os.path, threads: int, direct_io: bool, advise_random: bool):
		base = [
			"taskset",
			"0x1",
			exe,
			"--quiet",
			"oaoah",
			"-i", "{}".format(infile),
			"--dca", "{}".format(algo[:algo.find("-")]),
			"--query-file", "{}".format(queries),
			"-j", "{}".format(threads),
		]
		if direct_io:
			base.append("--direct-io")
		if advise_random:
			base.append("--advise-random")
		if args.checkpoint:
			base += ["--checkpoint-size", "{}".format(args.checkpoint)]
		return base

	@typechecked
	def run(self, num_queries: int, algo: str) -> None:
		io_opts : typing.List[bool] = []
		if args.dio:
			io_opts=[[True, False], [False, True], [False, False]]
		else:
			io_opts=[[False, True], [False, False]]

		for threads, [dio, adv_rnd] in itertools.product([1, 8], io_opts):
			dio: bool
			adv_rnd: bool
			if not os.path.exists(query_fn(num_queries)):
				logging.error("Query file {} does not exist".format(query_fn(num_queries)))
				sys.exit(-1)
			cmd_opts = {
					"exe":args.exe,
					"infile":infiles[algo],
					"algo":algo,
					"queries":query_fn(num_queries),
					"threads":threads,
					"direct_io": dio,
					"advise_random":adv_rnd
			}
			results.begin_run(algo, num_queries, threads, dio, adv_rnd)
			super().run(self.construct_cmdline(**cmd_opts))
			
			
class IMBench(Bench):

	@typechecked
	def construct_cmdline(self, exe: os.path, infile: os.path, algo: str, queries: os.path):
		base = [
			"taskset",
			"0x1",
			exe,
			"--quiet",
			"boost",
			"-i", "{}".format(infile),
			"--dca", "{}".format(algo[:algo.find("-")]),
			"--query-file", "{}".format(queries),
			"-j", "{}".format(1),
		]
		if args.checkpoint:
			base += ["--checkpoint-size", "{}".format(args.checkpoint)]
		return base

	@typechecked
	def run(self, num_queries: int, algo: str) -> None:
		if not os.path.exists(query_fn(num_queries)):
			logging.error("Query file {} does not exist".format(query_fn(num_queries)))
			sys.exit(-1)
		cmd_opts = {
				"exe":args.exe,
				"infile":infiles[algo],
				"algo":algo,
				"queries":query_fn(num_queries)
		}
		results.begin_run(algo, num_queries, 1, False, False)
		super().run(self.construct_cmdline(**cmd_opts))


# Do the actual benchmarks
for algo in ["dgc-oom", "chhl-oom", "ch-oom"]:
	if algo in args.bench:
		for num_queries in query_sizes:
			bench = OOMBench()
			bench.run(num_queries, algo)

for algo in ["dgc-im", "chhl-im", "ch-im"]:
	if algo in args.bench:
		for num_queries in query_sizes:
			bench = IMBench()
			bench.run(num_queries, algo)

if "sqloa" in args.bench:
	bench = DBBench(
		name=DBBench.dbname_from_filename(args.dgc_oom, True) if args.dbname is None else args.dbname,
		host=args.dbhost,
		port=args.dbport,
		user=args.dbuser,
		pw=args.dbpass,
		tablespace=args.dbspace
	)
	if args.dbimport:
		bench.import_data_oa(args.dgc_oom)
	
	for num_queries in query_sizes:
		bench.run_oa(num_queries)
	
	if args.dbdrop:
		bench.dropdb()

if "sql" in args.bench or "plsql" in args.bench:
	bench = DBBench(
		name=DBBench.dbname_from_filename(args.dgc_oom, False) if args.dbname is None else args.dbname,
		host=args.dbhost,
		port=args.dbport,
		user=args.dbuser,
		pw=args.dbpass,
		tablespace=args.dbspace
	)
	if args.dbimport:
		bench.import_data_bucketed(args.dgc_oom)
	
	if "sql" in args.bench:
		for num_queries in query_sizes:
			bench.run_sql(num_queries, False)
	if "plsql" in args.bench:
		for num_queries in query_sizes:
			bench.run_sql(num_queries, True)
	
	if args.dbdrop:
		bench.dropdb()
		
		
results.write_summary()

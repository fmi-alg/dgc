#include <iostream>
#include <vector>
#include <random>
#include <climits>
#include <fstream>
#include <queue>
#include <filesystem>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <sserialize/storage/UByteArrayAdapter.h>
#include <sserialize/Static/Array.h>
#include <sserialize/Static/Pair.h>
#include <sserialize/containers/OADHashTable.h>
#include <sserialize/algorithm/utilcontainerfuncs.h>
#include <sserialize/utility/assert.h>
#include <sserialize/stats/ProgressInfo.h>
#include <chrono>

#include <pqxx/pqxx>
#include <pqxx/version>

#ifndef NDEBUG
// 	#define VERBOSE
// 	#define VERBOSE_QUERY
//  #define GATHER_STATS
#endif

using clock_type = std::chrono::high_resolution_clock;

namespace c {
    constexpr int NO_ENTRY = -1;
}

struct DataBaseConfig {
	std::string host{"localhost"};
	std::string port{"5432"};
	std::string name;
	std::string user;
	std::string password;
	pqxx::connection connect() const {
		std::stringstream db_conn_params;
		db_conn_params << "host=" << host
						<< " port=" << port
						<< " user=" << user
						<< " password=" << password
						<< " dbname=" << name
						<< " connect_timeout=10";
		return pqxx::connection(db_conn_params.str());
	}
};

auto stream_to(pqxx::transaction_base & w, std::initializer_list<std::string_view> tblname, std::initializer_list<std::string_view> columns) {
	#if PQXX_VERSION_MAJOR < 7
	return pqxx::stream_to(w, std::string(*tblname.begin()), std::vector<std::string>(columns.begin(), columns.end()));
	#else
	return pqxx::stream_to::table(w, tblname, columns);
	#endif
}

auto write_values(auto & stream, auto... params) {
	#if PQXX_VERSION_MAJOR < 7
	stream << std::tuple<decltype(params)...>(params...);
	#else
	stream.write_values(std::forward<decltype(params)>(params)...);
	#endif
};

typedef int Rank;
typedef int Distance;
typedef std::pair<Rank, Distance> LabelElement;

static bool quiet = false;

std::ostream & operator<<(std::ostream & out, LabelElement const & el) {
	return out << "(" << el.first << ", " << el.second << ")";
}

namespace sserialize {
	template<>
	std::string nameOfType<::LabelElement>() { return "LabelElement"; }
}

struct StaticLabelElement: public std::pair<Rank, Distance> {
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
		ar & static_cast<std::pair<Rank, Distance> &>(*this);
	}
	StaticLabelElement() {}
	StaticLabelElement(sserialize::UByteArrayAdapter const & d) :
	std::pair<Rank, Distance>(
		d.get<Rank>(0),
		d.get<Distance>(sserialize::SerializationInfo<Rank>::length)
	)
	{}
};

class Bucket {
public:
	using Self = Bucket;
	using value_type = LabelElement;
	using const_iterator = sserialize::ReadOnlyAtStlIterator<Self*, value_type>;
	static constexpr sserialize::UByteArrayAdapter::SizeType EntrySize = sserialize::SerializationInfo<Rank>::length+sserialize::SerializationInfo<Distance>::length;
public:
	/// offset and size is in number of entries
	Bucket(std::size_t begin, std::size_t end, sserialize::UByteArrayAdapter const & d) :
	_off(begin),
	_size(end-begin),
	_d(d)
	{}
	~Bucket() {}
public:
	std::size_t size() const { return _size; }
	inline value_type at(std::size_t pos) const {
// 		if (pos >= size()) {
// 			throw std::out_of_range("");
// 		}
		return (*this)[pos];
	}
	inline value_type operator[](std::size_t pos) const {
		auto offset = (_off+pos)*EntrySize;
		return LabelElement(
			_d.get<Rank>(offset),
			_d.get<Distance>(offset+sserialize::SerializationInfo<Rank>::length)
		);
	}
	inline const_iterator begin() { return const_iterator(0, this); }
	inline const_iterator end() { return const_iterator(size(), this); }
private:
	std::size_t _off;
	std::size_t _size;
	sserialize::UByteArrayAdapter const & _d;
};

/// File format
// struct {
// 	uint64_t dataSize;
// 	UByteArrayAdapter _data;
// 	Array<uint32_t> _offsets;
// };
class LabelBuckets {
public:
	using Self = LabelBuckets;
	using value_type = Bucket;
	using const_iterator = sserialize::ReadOnlyAtStlIterator<Self*, value_type>;
	
public:
	LabelBuckets() {}
	LabelBuckets(LabelBuckets &&) = default;
	~LabelBuckets() {}
public:
	friend sserialize::UByteArrayAdapter & operator>>(sserialize::UByteArrayAdapter & src, LabelBuckets & dest);
public:
	inline std::size_t size() const { return _offsets.size()-1; }
	inline Bucket at(std::size_t pos) const {
		return (*this)[pos];
	}
	inline Bucket operator[](std::size_t pos) const {
		return Bucket(_offsets[pos], _offsets[pos+1], _data);
	}
public:
	inline const_iterator begin() { return const_iterator(0, this); }
	inline const_iterator end() { return const_iterator(size(), this); }
public:
	static void transform(std::vector<std::vector<LabelElement>> const & src, sserialize::UByteArrayAdapter & dest);
private:
	sserialize::UByteArrayAdapter _data;
	//Has an additional sentinel value
	sserialize::Static::Array<uint32_t> _offsets;
};

sserialize::UByteArrayAdapter & operator>>(sserialize::UByteArrayAdapter & src, LabelBuckets & dest) {
	#ifndef NDEBUG
	std::cerr << "Deserializing LabelBuckets at " << src.tellGetPtr() << std::endl;
	#endif
	uint64_t dataSize = src.getUint64();
	#ifndef NDEBUG
	std::cerr << "dataSize=" << dataSize << std::endl;
	#endif
	dest._data = sserialize::UByteArrayAdapter(src, src.tellGetPtr(), dataSize);
	SSERIALIZE_CHEAP_ASSERT_EQUAL(dataSize, dest._data.size());
	SSERIALIZE_CHEAP_ASSERT_EQUAL(0, dest._data.tellGetPtr());
	src.incGetPtr(dataSize);
	#ifndef NDEBUG
	std::cerr << "offset index starts at " << src.tellGetPtr() << std::endl;
	#endif
	src >> dest._offsets;
	#ifndef NDEBUG
	std::cerr << "Deserializing LabelBucket finished at " << src.tellGetPtr() << std::endl;
	#endif
	return src;
}

void LabelBuckets::transform(const std::vector<std::vector<LabelElement>> & src, sserialize::UByteArrayAdapter& dest)
{
	std::cerr << "Serializing LabelBuckets at " << dest.tellPutPtr() << std::endl;
	std::size_t dataSizePtr = dest.tellPutPtr();
	dest.putUint64(0);
	std::size_t dataBeginPtr = dest.tellPutPtr();
	std::vector<uint32_t> offsets;
	offsets.reserve(src.size());
	std::size_t count = 0;
	for(auto const & bucket : src) {
		offsets.push_back(count);
		for(auto const & label : bucket) {
			dest.put<Rank>(label.first);
			dest.put<Distance>(label.second);
		}
		count += bucket.size();
	}
	//The sentinel value
	offsets.push_back(count);
	if (count > std::numeric_limits<uint32_t>::max()) {
		throw std::out_of_range("Too many labels");
	}
	auto dataSize = dest.tellPutPtr()-dataBeginPtr;
	SSERIALIZE_CHEAP_ASSERT_EQUAL((count-1)*Bucket::EntrySize, dataSize);
	dest.putUint64(dataSizePtr, dataSize);
	#ifndef NDEBUG
	std::cerr << "dataSize=" << dataSize << std::endl;
	std::cerr << "offset index starts at " << dest.tellPutPtr() << std::endl;
	#endif
	dest << offsets;
	#ifndef NDEBUG
	std::cerr << "Serializing LabelBucket finished at " << dest.tellPutPtr() << std::endl;
	#endif
}

/// Labels table: (position serial, tgt int not null, dst int not null)
/// Offset table: (position serial, begin int not null, end int not null)
class DBLabels {
public:
	using Self = DBLabels;
	enum Type {
		Forward, Backward
	};
	struct QueryResult {
		QueryResult(pqxx::result && result) : result(std::move(result)), iter(this->result.iter<Rank, Distance>()) {}
		QueryResult(QueryResult const &) = delete;
		QueryResult(QueryResult && other) : result(std::move(result)), iter(this->result.iter<Rank, Distance>()) {}
		pqxx::result result;
		pqxx::internal::result_iteration<Rank, Distance> iter;
		auto begin() const { return iter.begin(); }
		auto end() const { return iter.end(); }
	};
	using value_type = std::vector<QueryResult>;
	using const_iterator = sserialize::ReadOnlyAtStlIterator<Self*, value_type>;
	
public:
	DBLabels(DataBaseConfig const & cfg, Type type, std::size_t num_labels = 0);
	DBLabels(DBLabels &&) = default;
	~DBLabels() {
#ifdef GATHER_STATS
		std::cout << "Execution times for ";
		if (_type == Type::Forward) {
			std::cout << "forward";
		}
		else {
			std::cout << "backward";
		}
		std::cout << " labels [us]:" << std::endl;
		std::chrono::microseconds total{0};
		for(auto [nodeid, time] : _query_times) {
			std::cout << nodeid << ": " << time.count() << std::endl;
			total += time;
		}
		std::cout << "Total: " << total.count() << std::endl;
#endif
	}
public:
	std::size_t size() const { return _size; }
	auto at(std::size_t pos) const {
#ifdef GATHER_STATS
		auto start_time = clock_type::now();
#endif
		pqxx::nontransaction ta(_con);
		value_type result;
#if 0
		int ofs_begin, ofs_end;
		{
			auto r = ta.exec_prepared1("get_offsets", pos);
			ofs_begin = r[0].as<int>();
			ofs_end = r[1].as<int>();
		}
		{
			auto r = ta.exec_prepared("get_label_from_offsets", ofs_begin, ofs_end);
			for(auto row : r) {
				result.emplace_back(row[0].as<int>(), row[1].as<int>());
			}
		}
#elif 0
		auto r = ta.exec_prepared("get_label_from_nodeid", pos);
		for(auto row : r) {
			result.emplace_back(row[0].as<int>(), row[1].as<int>());
		}
#elif 0
		auto r = ta.exec("select * from " + std::string(qfname()) + "(" + std::to_string(pos) + ");");
		for(auto row : r) {
			result.emplace_back(row[0].as<int>(), row[1].as<int>());
		}
#elif 1
		auto r = ta.exec("select * from " + std::string(qfname()) + "(" + std::to_string(pos) + ");");
#ifdef GATHER_STATS
		_query_times[pos] = std::chrono::duration_cast<std::chrono::microseconds>(clock_type::now()-start_time);
#endif
		return QueryResult(std::move(r));
#elif 1
//does not work since COPY does not allow EXECUTE
// 		return ta.stream<Rank, Distance>("EXECUTE get_label_from_nodeid(" + std::to_string(pos) + ")");
#endif
// 		return result;
	}
	inline auto operator[](std::size_t pos) const { return at(pos); }
public:
	inline const_iterator begin() { return const_iterator(0, this); }
	inline const_iterator end() { return const_iterator(size(), this); }
public:
	std::string_view qfname() const {
		if (_type == Forward) {
			return "fw_get_label";
		}
		else {
			return "bw_get_label";
		}
	}
	std::string_view tablename() const {
		if (_type == Forward) {
			return "fwlabels";
		}
		else {
			return "bwlabels";
		}
	}
	std::string_view offsetsname() const {
		if (_type == Forward) {
			return "fwoffsets";
		}
		else {
			return "bwoffsets";
		}
	}
public:
	mutable pqxx::connection _con;
	Type _type;
	std::size_t _size;
#ifdef GATHER_STATS
	mutable std::unordered_map<std::size_t, std::chrono::microseconds> _query_times;
#endif
};

DBLabels::DBLabels(DataBaseConfig const & cfg, Type t, std::size_t num_labels) :
_con(cfg.connect()),
_type(t)
{
	//prepare full query
	{
		std::stringstream ss;
		ss <<
	"select \
		l.tgt, l.dst \
	from\n"
	<< tablename() << " as l, \
		( \
			select \
				ofs.obegin, \
				ofs.oend\
			from\n"
	<< offsetsname() << " as ofs \
			where \
				ofs.pos = $1 \
			limit 1 \
		) as ofs \
	where \
		ofs.obegin <= l.pos and l.pos < ofs.oend \
	;";
// 		std::cout << ss.str() << std::endl;
		_con.prepare("get_label_from_nodeid", ss.str());
	}
	//prepare offset select
	{
		std::stringstream ss;
		ss <<
		"select\n"
		"\tofs.obegin as b,\n"
		"\tofs.oend as e\n"
		"from\n"
		<< offsetsname() << " as ofs\n"
		"where\n"
		"\tofs.pos = $1\n"
		"limit 1;";
// 		std::cout << ss.str() << std::endl;
		_con.prepare("get_offsets", ss.str());
	}
	//prepare label retrieval from offsets
	{
		std::stringstream ss;
		ss <<
		"select\n"
		"\tl.tgt, l.dst\n"
		"from\n"
		"\t" << tablename() << " as l\n"
		"where\n"
		"\t$1 <= l.pos and l.pos < $2;";
		_con.prepare("get_label_from_offsets", ss.str());
	}
	if (!num_labels) { //get number of nodes
		try {
			pqxx::nontransaction t(_con);
			_size = t.query_value<std::size_t>("select stats.node_count from stats limit 1;");
		}
		catch(std::exception const &) {
#ifdef VERBOSE
			std::cout << "Counting number of labels" << std::endl;
#endif
			pqxx::nontransaction t(_con);
			std::string stmt = "select count(*) from " + std::string(offsetsname()) + ";";
			_size = t.query_value<std::size_t>(stmt);
		}
	}
	//Create query function for single roundtrip fetch
	{
		pqxx::nontransaction t(_con);
		std::stringstream ss;
		ss << "CREATE OR REPLACE FUNCTION " << qfname() << "(IN nodeid int)\n";
		ss <<
R"...(
		RETURNS table (
			nid int,
			dst int
		)
		AS
		$$
	DECLARE
		ofs_begin bigint := 0;
		ofs_end bigint := 0;
	begin
		select
			ofs.obegin, ofs.oend
		into
			ofs_begin, ofs_end
		from
			)..."
		<< offsetsname() << " as ofs\n"
R"...(		where
			ofs.pos = nodeid
		limit 1;

		return query
			select
				fw.tgt, fw.dst
			from
				)..."
			<< tablename() << " as fw"
R"...(
			where
				ofs_begin <= fw.pos and fw.pos < ofs_end
			;
	END
	$$ LANGUAGE plpgsql;)...";
		t.exec0(ss.str());
	}
}

namespace sserialize {

template<>
struct SerializationInfo<LabelElement> {
	static const bool is_fixed_length = SerializationInfo<Rank>::is_fixed_length && SerializationInfo<Distance>::is_fixed_length;
	static const OffsetType length = SerializationInfo<Rank>::length + SerializationInfo<Distance>::length;
	static const OffsetType max_length = SerializationInfo<Rank>::max_length + SerializationInfo<Distance>::max_length;
	static const OffsetType min_length = SerializationInfo<Rank>::min_length + SerializationInfo<Distance>::min_length;
	static OffsetType sizeInBytes(const LabelElement & value) {
		return length;
	}
};

template<>
struct SerializationInfo<StaticLabelElement>: SerializationInfo<LabelElement> {};

} //end namespace sserialize

sserialize::UByteArrayAdapter & operator<<(sserialize::UByteArrayAdapter & dest, LabelElement const & src) {
	return dest << src.first << src.second;
}

sserialize::UByteArrayAdapter & operator>>(sserialize::UByteArrayAdapter & src, LabelElement & dest) {
	return src >> dest.first >> dest.second;
}

sserialize::UByteArrayAdapter & operator>>(sserialize::UByteArrayAdapter & src, StaticLabelElement & dest) {
	return src >> dest.first >> dest.second;
}

class DistanceInfoVector {
public:
	using size_type = uint32_t;
	using key_type = uint32_t;
	using mapped_type = int;
	static constexpr mapped_type ne = c::NO_ENTRY;
public:
	class Iterator {
	public:
		using value_type = std::pair<key_type, mapped_type>;
	public:
		Iterator(DistanceInfoVector const * d, std::vector<key_type>::const_iterator it) :
		_d(d),
		_it(it)
		{}
		Iterator(Iterator const &) = default;
		Iterator(Iterator &&) = default;
		~Iterator() {}
		Iterator & operator=(Iterator const &) = default;
		Iterator & operator=(Iterator &&) = default;
	public:
		inline value_type operator*() const {
			return value_type(*_it, _d->_d[*_it]);
		}
		inline Iterator & operator++() {
			++_it;
			return *this;
		}
	public:
		friend inline bool operator!=(Iterator const & a, Iterator const & b) { return a._it != b._it; }
		friend inline bool operator==(Iterator const & a, Iterator const & b) { return a._it == b._it; }
	private:
		DistanceInfoVector const * _d;
		std::vector<key_type>::const_iterator _it;
	};
	using const_iterator = Iterator;
public:
	DistanceInfoVector() {}
	~DistanceInfoVector() {}
public:
	inline std::size_t size() const { return _ocpied.size(); }
	inline void reserve(size_type size) { _d.resize(size, ne); }
	void clear() {
		for(auto x : _ocpied) {
			_d[x] = ne;
		}
		_ocpied.clear();
	}
	inline void emplace(key_type key, mapped_type value) {
		if (_d.at(key) != ne) {
			throw std::runtime_error("Key already inserted");
		}
		_ocpied.push_back(key);
		_d.at(key) = value;
	}
	inline bool count(key_type key) {
		return _d.at(key) != ne;
	}
	inline const_iterator begin() const {
		return const_iterator(this, _ocpied.begin());
	}
	inline const_iterator end() const {
		return const_iterator(this, _ocpied.end());
	}
private:
	inline void check_access(key_type key) const {
		#ifndef NDEBUG
		if (_d.at(key) == ne) {
			throw std::runtime_error("Value of key is not valid");
		}
		#endif
	}
public: //The following functions don't track insertions!
	inline mapped_type const & operator[](key_type key) const {
		check_access(key);
		return _d[key];
	}
	inline mapped_type & operator[](key_type key) {
		check_access(key);
		return _d[key];
	}
	inline mapped_type const & at(key_type key) const {
		check_access(key);
		return _d.at(key);
	}
	inline mapped_type & at(key_type key) {
		check_access(key);
		return _d.at(key);
	}
private:
	std::vector<mapped_type> _d;
	std::vector<key_type> _ocpied;
};

enum class DistanceComputationAlgo {
	HopGroups,
	CHHL,
	CH
};

std::mt19937 my_rng;

template<
	typename TForwardBuckets,
	typename TBackwardBuckets,
	typename TGroupBounds,
	typename TDistanceInfo = DistanceInfoVector
>
class DGCLabeling
{
public:
	using ForwardBuckets = TForwardBuckets;
	using BackwardBuckets = TBackwardBuckets;
	using GroupBounds = TGroupBounds;
	using DistanceInfo = TDistanceInfo;
public:
	DGCLabeling() {}
	DGCLabeling(DGCLabeling &&) = default;
	DGCLabeling(ForwardBuckets && fwb, BackwardBuckets && bwb, GroupBounds && gb) :
	_fwdLabels(std::move(fwb)),
	_bwdLabels(std::move(bwb)),
	_groupBounds(std::move(gb))
	{
		_numNodes = _fwdLabels.size();
        _numHops = _groupBounds.size();
	}
	void stats() {
		std::cout << "nodes: " << _numNodes << '\n';
		std::cout << "hops:" << _numHops << '\n';
	}
	GroupBounds const & groubBounds() const { return _groupBounds; }
	std::size_t numNodes() const { return _numNodes; }
	
	/// @return .second has the computation time for the first .first queries
	/// @param checkpoint_size how many queries should be computed before we take a measurement
	std::vector<std::pair<std::size_t, std::chrono::microseconds>> speedTest(std::vector<std::pair<int, int>> const & queries, DistanceComputationAlgo dca, std::size_t threadCount, std::size_t checkpoint_size) {
		if (!checkpoint_size) {
			checkpoint_size = queries.size();
		}
		long sumDistances = 0;
		std::vector<std::pair<std::size_t, std::chrono::microseconds>> cummulative_times;
		cummulative_times.reserve(queries.size()/checkpoint_size + (queries.size()%checkpoint_size!=0) + 1);
		clock_type::time_point start_time;
		if (threadCount > 1) {
			std::vector<std::thread> threads;
			threads.reserve(threadCount);
			std::atomic<std::size_t> i{0};
			std::atomic<std::size_t> finished{0};
			std::atomic<long> a_sumDistances{0};
			std::atomic<std::size_t> setup_complete{0};
			std::atomic<bool> wait_for_go{true};
			std::mutex lck;
			for(std::size_t j(0); j < threadCount; ++j) {
				threads.emplace_back([&,this](){
					DistanceComputation dc(this);
					long my_sumDistances = 0;
					setup_complete.fetch_add(1, std::memory_order_seq_cst);
					while(wait_for_go.load(std::memory_order_seq_cst)) { //do some very nasty spinning
						;
					}
					while (true) {
						auto qi = i.fetch_add(1, std::memory_order_seq_cst);
						if (qi >= queries.size()) {
							break;
						}
						switch(dca) {
							case DistanceComputationAlgo::HopGroups:
								my_sumDistances += dc.getDistance(queries[qi].first, queries[qi].second);
								break;
							case DistanceComputationAlgo::CHHL:
								my_sumDistances += dc.template getDistanceCHHL<false>(queries[qi].first, queries[qi].second);
								break;
							case DistanceComputationAlgo::CH:
								my_sumDistances += dc.template getDistanceCHHL<true>(queries[qi].first, queries[qi].second);
								break;
						}
						auto finish_count = finished.fetch_add(1, std::memory_order_seq_cst)+1;
						if (finish_count && finish_count%checkpoint_size == 0) {
							auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(clock_type::now()-start_time);
							std::lock_guard<std::mutex> l(lck);
							cummulative_times.emplace_back(finish_count+1, elapsed_time);
						}
					}
					a_sumDistances.fetch_add(my_sumDistances, std::memory_order_relaxed);
				});
			}
			while(setup_complete.load(std::memory_order_seq_cst) != threadCount) { //spin until all threads are ready for action
				;
			}
			start_time = clock_type::now();
			wait_for_go.store(false, std::memory_order_seq_cst);
			for(auto & x : threads) {
				x.join();
			}
			sumDistances = a_sumDistances;
		}
		else {
			DistanceComputation dc(this);
			start_time = clock_type::now();
			for (std::size_t i(0), s(queries.size()); i < s;) {
				for(std::size_t j(0); j < checkpoint_size && i < s; ++j, ++i) {
					switch(dca) {
					case DistanceComputationAlgo::HopGroups:
						sumDistances += dc.getDistance(queries[i].first, queries[i].second);
						break;
					case DistanceComputationAlgo::CHHL:
						sumDistances += dc.template getDistanceCHHL<false>(queries[i].first, queries[i].second);
						break;
					case DistanceComputationAlgo::CH:
						sumDistances += dc.template getDistanceCHHL<true>(queries[i].first, queries[i].second);
						break;
					}
				}
				cummulative_times.emplace_back(i, std::chrono::duration_cast<std::chrono::microseconds>(clock_type::now()-start_time));
			}
		}
		std::cout << "sum of distances: " << sumDistances << std::endl;
        return cummulative_times;
    }
    std::vector<std::pair<int, int>> gen_queries(std::size_t count) const {
        std::vector<std::pair<int, int>> queries;
		auto d = std::uniform_int_distribution<int>(0, _numNodes-1);
		for(std::size_t i(0); i < count; ++i) {
			queries.emplace_back(d(my_rng), d(my_rng));
		}
		return queries;
	}
    std::vector<std::pair<std::size_t, std::chrono::microseconds>> speedTest(int numRuns, DistanceComputationAlgo dca, std::size_t threadCount, std::size_t checkpoint_size) {
        return speedTest(gen_queries(numRuns), dca, threadCount, checkpoint_size);
	}
	///@return (distance, time in nsecs)
	auto getTimedDistance(int source, int target) {
		DistanceComputation dc(this);
		auto begin = std::chrono::high_resolution_clock::now();
		auto dist = dc.getDistance(source, target);
		auto end = std::chrono::high_resolution_clock::now();
		auto elapsed = end-begin;
		return std::make_pair(dist, std::chrono::duration_cast<std::chrono::microseconds>(elapsed));
	}
	
	void export_postgresql_stats(pqxx::nontransaction & w) {
		std::string stmt = "DROP TABLE IF EXISTS stats; CREATE UNLOGGED TABLE stats (group_count bigint not null, node_count bigint not null);";
		std::cout << stmt << std::endl;
		w.exec(stmt);
		{
			auto stream = stream_to(w, {"stats"}, {"group_count", "node_count"});
			write_values(stream, _numHops, _numNodes);
			stream.complete();
		}
		stmt = "DROP TABLE IF EXISTS bounds; CREATE UNLOGGED TABLE bounds (grp bigint not null, bound bigint not null);";
		std::cout << stmt << std::endl;
		w.exec(stmt);
		{
			auto stream = stream_to(w, {"bounds"}, {"grp", "bound"});
			for(std::size_t i(0); i < _groupBounds.size(); ++i) {
				write_values(stream, i, _groupBounds.at(i));
			}
			stream.complete();
		}
	}
	
	//Store data as offset array in postgresdb
	void export_postgresql_oa(pqxx::connection & conn) {
		pqxx::nontransaction w(conn);
		
		auto export_labels = [&](std::string const & prefix, auto const & buckets) {
			sserialize::ProgressInfo pinfo;
			std::string stmt;
			std::string labeltbl = prefix + "labels";
			std::string offsetbl = prefix + "offsets";
			
			stmt = "DROP TABLE IF EXISTS " + labeltbl + "; CREATE UNLOGGED TABLE " + labeltbl + " (pos bigint not null, tgt int not null, dst int not null);";
			std::cout << stmt << std::endl;
			w.exec(stmt);
			std::vector<std::pair<std::size_t, std::size_t>> offsets;
			offsets.reserve(buckets.size());
			{
				auto labels_stream = stream_to(w, {labeltbl}, {"pos", "tgt", "dst"});
				std::size_t offset = 0;
				pinfo.begin(buckets.size(), labeltbl);
				for(std::size_t src(0), s(buckets.size()); src < s; ++src) {
					std::size_t offset_begin = offset;
					for(LabelElement le : buckets[src]) {
						write_values(labels_stream, offset, le.first, le.second);
						++offset;
					}
					offsets.emplace_back(offset_begin, offset);
					pinfo(src);
				}
				labels_stream.complete();
				pinfo.end();
			}
			stmt = "create unique index " + labeltbl + "_idx on " + labeltbl + " using btree(pos asc nulls last);";
			std::cout << stmt << std::endl;
			w.exec(stmt);
			
			//end is one passed the end
			stmt = "DROP TABLE IF EXISTS "+ offsetbl + "; CREATE UNLOGGED TABLE " + offsetbl +" (pos int not null, obegin bigint not null, oend bigint not null);";
			std::cout << stmt << std::endl;
			w.exec(stmt);
			{
				auto offsets_stream = stream_to(w, {offsetbl}, {"pos", "obegin", "oend"});
				std::size_t count = 0;
				pinfo.begin(offsets.size(), offsetbl);
				for(auto [b, e] : offsets) {
					write_values(offsets_stream, count, b, e);
					pinfo(count++);
				}
				offsets_stream.complete();
				pinfo.end();
			}
			stmt = "create unique index " + offsetbl + "_idx on " + offsetbl + " using btree(pos asc nulls last);";
			std::cout << stmt << std::endl;
			w.exec(stmt);
		};
		export_labels("fw", _fwdLabels);
		export_labels("bw", _bwdLabels);
		std::string stmt = "DROP TABLE IF EXISTS grp; CREATE UNLOGGED TABLE grp (pos bigint not null, grp int not null);";
		std::cout << stmt << std::endl;
		w.exec(stmt);
		{
			auto group_stream = stream_to(w, {"grp"}, {"pos", "grp"});
			for(std::size_t src(0), s(_fwdLabels.size()); src < s; ++src) {
				write_values(group_stream, src, _getGroupOfNode(src));
			}
			group_stream.complete();
		}
		stmt = "create unique index grp_idx on grp using btree(pos asc nulls last);";
		std::cout << stmt << std::endl;
		w.exec(stmt);
		export_postgresql_stats(w);
		w.commit();
	}
	//One table per group, group info in extra table, entries not null
	//index only on src since the other columns are not needed
	//ascending and nulls last
	void export_postgresql(pqxx::connection & conn) {
		pqxx::nontransaction w(conn);
		
		sserialize::ProgressInfo pinfo;
		for(std::size_t grp(0); grp < _numHops; ++grp) {
			std::string fwtbln = "fwdlabels" + std::to_string(grp);
			std::string bwtbln = "bwdlabels" + std::to_string(grp);
			w.exec("DROP TABLE IF EXISTS " + fwtbln + "; CREATE UNLOGGED TABLE " + fwtbln + "(src int not null, tgt int not null, dst int not null)");
			w.exec("DROP TABLE IF EXISTS " + bwtbln + "; CREATE UNLOGGED TABLE " + bwtbln + "(src int not null, tgt int not null, dst int not null)");
			{
				auto stream = stream_to(w, {fwtbln}, {"src", "tgt", "dst"});
				pinfo.begin(_fwdLabels.size(), "Forward labels group " + std::to_string(grp));
				for(std::size_t src(0), s(_fwdLabels.size()); src < s; ++src) {
					if (_getGroupOfNode(src) == grp) {
						for(LabelElement le : _fwdLabels[src]) {
							write_values(stream, src, le.first, le.second);
						}
					}
					pinfo(src);
				}
				stream.complete();
				pinfo.end();
			}
			{
				auto stream = stream_to(w, {bwtbln}, {"src", "tgt", "dst"});
				pinfo.begin(_bwdLabels.size(), "Backward labels group " + std::to_string(grp));
				for(std::size_t src(0), s(_bwdLabels.size()); src < s; ++src) {
					if (_getGroupOfNode(src) == grp) {
						for(LabelElement le : _bwdLabels[src]) {
							write_values(stream, src, le.first, le.second);
						}
					}
					pinfo(src);
				}
				stream.complete();
				pinfo.end();
			}
		}
		{
			assert(_bwdLabels.size() == _fwdLabels.size());
			w.exec("DROP TABLE IF EXISTS grp; CREATE UNLOGGED TABLE grp (node serial, grp int)");
			auto stream = stream_to(w, {"grp"}, {"grp"});
			pinfo.begin(_bwdLabels.size(), "Group info");
			for(std::size_t src(0), s(_bwdLabels.size()); src < s; ++src) {
				write_values(stream, _getGroupOfNode(src));
				pinfo(src);
			}
			stream.complete();
			pinfo.end();
		}
		{ //create index
			if (!quiet) {
				std::cout << "Creating indexes" << std::endl;
			}
			pinfo.begin(_numHops, "Indexes");
			for(std::size_t grp(0); grp < _numHops; ++grp) {
				std::string fwtbln = "fwdlabels" + std::to_string(grp);
				std::string bwtbln = "bwdlabels" + std::to_string(grp);
				w.exec("CREATE INDEX " + fwtbln + "_src_idx ON " + fwtbln + " USING btree (src ASC)");
				w.exec("CREATE INDEX " + bwtbln + "_src_idx ON " + bwtbln + " USING btree (src ASC)");
				pinfo(grp);
			}
			w.exec("CREATE INDEX grp_node_idx ON grp USING btree (node ASC)");
			pinfo.end();
		}
		export_postgresql_stats(w);
		w.commit();
	}
	
	struct DistanceComputation {
		DistanceComputation(DGCLabeling * parent) : parent(parent) {
			if constexpr (std::is_same_v<DistanceInfoVector, DistanceInfo>) {
				_fwdDistances.reserve(parent->_numNodes);
				_bwdDistances.reserve(parent->_numNodes);
			}
			else {
				if constexpr (!std::is_same_v<std::map<uint32_t, int>, DistanceInfo>) { //Let's reserve at least 16 MBytes
					std::size_t reserve_count = 16*1024*1024/sizeof(std::pair<Rank, Distance>);
					_fwdDistances.reserve(reserve_count);
					_bwdDistances.reserve(reserve_count);
				}
			}
		}
		
		int getDistance(int source, int target)
		{
			if (source == target)
			{
				return 0;
			}
#ifndef NDEBUG
			if (source < 0 || parent->_numNodes <= source)
			{
				throw std::runtime_error("Source node id is out of range.");
			}
			if (target < 0 || parent->_numNodes <= target)
			{
				throw std::runtime_error("Target node id is out of range.");
			}
#endif
	#ifdef VERBOSE_QUERY
		{
			auto bucket = parent->_fwdLabels[source];
			sserialize::print<','>(std::cout, bucket.begin(), bucket.end());
		}
	#endif
			std::vector<std::vector<int>> fwd_buckets(parent->_numHops);
			std::vector<std::vector<int>> bwd_buckets(parent->_numHops);
			int source_group = parent->_getGroupOfNode(source);
			_fwdDistances.emplace(source, 0);
			for (auto [nodeid, dist] : parent->_fwdLabels[source])
			{
				int new_distance = dist;
				if (!_fwdDistances.count(nodeid))
				{
					_fwdDistances.emplace(nodeid, new_distance);
					int next_group = parent->_getGroupOfNode(nodeid);
					if (source_group < next_group)
					{
						fwd_buckets[next_group].push_back(nodeid);
					}
				}
				else if (new_distance < _fwdDistances[nodeid])
				{
					_fwdDistances[nodeid] = new_distance;
				}
			}
			for (int group = source_group + 1; group < parent->_numHops; group++)
			{
				for (int rank : fwd_buckets[group])
				{
					int distance = _fwdDistances[rank];
					for (auto [nodeid,  dist] : parent->_fwdLabels[rank])
					{
						int new_distance = dist + distance;
						if (!_fwdDistances.count(nodeid))
						{
							_fwdDistances.emplace(nodeid, new_distance);
							int next_group = parent->_getGroupOfNode(nodeid);
							if (group < next_group)
							{
								fwd_buckets[next_group].push_back(nodeid);
							}
						}
						else if (new_distance < _fwdDistances[nodeid])
						{
							_fwdDistances[nodeid] = new_distance;
						}
					}
				}
			}
			int target_group = parent->_getGroupOfNode(target);
			_bwdDistances.emplace(target, 0);
			for (auto [nodeid, dist] : parent->_bwdLabels[target])
			{
				int new_distance = dist;
				if (!_bwdDistances.count(nodeid))
				{
					_bwdDistances.emplace(nodeid, new_distance);
					int next_group = parent->_getGroupOfNode(nodeid);
					if (target_group < next_group)
					{
						bwd_buckets[next_group].push_back(nodeid);
					}
				}
				else if (new_distance < _bwdDistances[nodeid])
				{
					_bwdDistances[nodeid] = new_distance;
				}
			}
			for (int group = target_group + 1; group < parent->_numHops; group++)
			{
				for (int rank : bwd_buckets[group])
				{
					int distance = _bwdDistances[rank];
					for (auto [nodeid, dist] : parent->_bwdLabels[rank])
					{
						int new_distance = dist + distance;
						if (!_bwdDistances.count(nodeid))
						{
							_bwdDistances.emplace(nodeid, new_distance);
							int next_group = parent->_getGroupOfNode(nodeid);
							if (group < next_group)
							{
								bwd_buckets[next_group].push_back(nodeid);
							}
						}
						else if (new_distance < _bwdDistances[nodeid])
						{
							_bwdDistances[nodeid] = new_distance;
						}
					}
				}
			}
			#ifdef VERBOSE_QUERY
			std::cout << "Forward buckets:\n";
			for(int bucket(0); bucket < parent->_numHops; ++bucket) {
				std::cout << bucket << ": ";
				sserialize::print(std::cout, fwd_buckets[bucket]);
				std::cout << std::endl;
			}
			std::cout << "Forward distances (" << _fwdDistances.size() << "):";
			#endif
			int best_dist = c::NO_ENTRY;
			for (auto [rank, distance] : _fwdDistances)
			{
	#ifdef VERBOSE_QUERY
				std::cout << "(" << rank << ", " << distance << ", " << parent->_getGroupOfNode(rank) << ")\n";
	#endif
				if (_bwdDistances.count(rank))
				{
					int new_distance = _bwdDistances[rank] + distance;
					if (best_dist == c::NO_ENTRY || new_distance < best_dist)
					{
						best_dist = new_distance;
					}
				}
			}
			_fwdDistances.clear();
			_bwdDistances.clear();
			return best_dist;
		}
		template<bool CHLike>
		int getDistanceCHHL(int source, int target)
		{
			if (!CHLike && parent->_numHops > 2) {
				throw std::runtime_error("Too many groups");
			}
			if (source == target)
			{
				return 0;
			}
#ifndef NDEBUG
			if (source < 0 || parent->_numNodes <= source)
			{
				throw std::runtime_error("Source node id is out of range.");
			}
			if (target < 0 || parent->_numNodes <= target)
			{
				throw std::runtime_error("Target node id is out of range.");
			}
#endif
			std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pQ;
			pQ.push(std::make_pair(0, source));
			_fwdDistances.emplace(source, 0);
			while (!pQ.empty())
			{
				std::pair<int, int> next_element = pQ.top();
				int rank = next_element.second;
				int distance = next_element.first;
				pQ.pop();
				if (_fwdDistances[rank] == distance)
				{
					bool should_be_stalled = false;
					for (auto [nodeid, dst] : parent->_bwdLabels[rank])
					{
						if (_fwdDistances.count(nodeid))
						{
							int new_dist = _fwdDistances[nodeid] + dst;
							if (new_dist < distance)
							{
								should_be_stalled = true;
								break;
							}
						}
					}
					if (!should_be_stalled)
					{
	                    int group;
						if (!CHLike) {
							group = parent->_getGroupOfNode(rank);
						}
						for (auto [nodeid, dist] : parent->_fwdLabels[rank])
						{
							int new_distance = dist + distance;
							if (!_fwdDistances.count(nodeid))
							{
								_fwdDistances.emplace(nodeid, new_distance);
	                            if (CHLike || group == 0)
								{
									pQ.push(std::make_pair(new_distance, nodeid));
								}
							}
							else if (new_distance < _fwdDistances[nodeid])
							{
								_fwdDistances[nodeid] = new_distance;
	                            if (CHLike || group == 0)
								{
									pQ.push(std::make_pair(new_distance, nodeid));
								}
							}
						}
					}
				}
			}
			pQ.push(std::make_pair(0, target));
			_bwdDistances.emplace(target, 0);
			while (!pQ.empty())
			{
				std::pair<int, int> next_element = pQ.top();
				int rank = next_element.second;
				int distance = next_element.first;
				pQ.pop();
				if (_bwdDistances[rank] == distance)
				{
					bool should_be_stalled = false;
					for (auto [nodeid, dist] : parent->_fwdLabels[rank])
					{
						if (_bwdDistances.count(nodeid))
						{
							int new_dist = _bwdDistances[nodeid] + dist;
							if (new_dist < distance)
							{
								should_be_stalled = true;
								break;
							}
						}
					}
					if (!should_be_stalled)
					{
	                    int group;
						if (!CHLike) {
							group = parent->_getGroupOfNode(rank);
						}
						for (auto [nodeid, dist] : parent->_bwdLabels[rank])
						{
							int new_distance = dist + distance;
							if (!_bwdDistances.count(nodeid))
							{
								_bwdDistances.emplace(nodeid, new_distance);
	                            if (CHLike || group == 0)
								{
									pQ.push(std::make_pair(new_distance, nodeid));
								}
							}
							else if (new_distance < _bwdDistances[nodeid])
							{
								_bwdDistances[nodeid] = new_distance;
	                            if (CHLike || group == 0)
								{
									pQ.push(std::make_pair(new_distance, nodeid));
								}
							}
						}
					}
				}
			}
			int best_dist = c::NO_ENTRY;
			for (auto [rank, distance] : _fwdDistances)
			{
				if (_bwdDistances.count(rank))
				{
					int new_distance = _bwdDistances[rank] + distance;
					if (best_dist == c::NO_ENTRY || new_distance < best_dist)
					{
						best_dist = new_distance;
					}
				}
			}
			_fwdDistances.clear();
			_bwdDistances.clear();

			return best_dist;
		}
		DGCLabeling * parent;
		DistanceInfo _fwdDistances;
		DistanceInfo _bwdDistances;
	};
    

protected:
	//Stuff that is serialized
    ForwardBuckets _fwdLabels;
    BackwardBuckets _bwdLabels;
    GroupBounds _groupBounds;
	
	//dynamic stuff
    int _numNodes{0};
    int _numHops{0};

    virtual int _getGroupOfNode(int rank)
    {
        int group = 0;
#ifdef NDEBUG
        while (rank >= _groupBounds[group])
#else
        while (rank >= _groupBounds.at(group))
#endif
        {
            group++;
        }
        return group;
    }
};

class DGCLabeling_Boost: public
	DGCLabeling<
		std::vector<std::vector<LabelElement>>,
		std::vector<std::vector<LabelElement>>,
		std::vector<int>
	>
{
public:
	using Self = DGCLabeling_Boost;
	using MyBaseClass = 	DGCLabeling<
		std::vector<std::vector<LabelElement>>,
		std::vector<std::vector<LabelElement>>,
		std::vector<int>
	>;
	using ForwardBuckets = MyBaseClass::ForwardBuckets;
	using BackwardBuckets = MyBaseClass::BackwardBuckets;
	using GroupBounds = MyBaseClass::GroupBounds;
public:
	struct SerializationWrapper {
		ForwardBuckets & fwb;
		BackwardBuckets & bwb;
		GroupBounds & gb;
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int version)
		{
			ar & fwb;
			ar & bwb;
			ar & gb;
		}
	};
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
		SerializationWrapper sw{
			.fwb = _fwdLabels,
			.bwb = _bwdLabels,
			.gb = _groupBounds
		};
		ar & sw;
    }

    static Self fromBinary(std::string file)
    {
        std::ifstream myHLfile(file, std::ios_base::in | std::ios_base::binary);
        boost::archive::binary_iarchive ia(myHLfile);
        ForwardBuckets fwb;
		BackwardBuckets bwb;
		GroupBounds gb;
		SerializationWrapper sw{
			.fwb = fwb,
			.bwb = bwb,
			.gb = gb
		};
		ia & sw;
		return Self(std::move(fwb), std::move(bwb), std::move(gb));
    }
    
    void sserialize(sserialize::UByteArrayAdapter & dest) const {
		dest << _fwdLabels << _bwdLabels << _groupBounds;
	}
	
	void sserialize_oa(sserialize::UByteArrayAdapter & dest) const {
		LabelBuckets::transform(_fwdLabels, dest);
		LabelBuckets::transform(_bwdLabels, dest);
		dest << _groupBounds;
	}
	DGCLabeling_Boost() {}
	DGCLabeling_Boost(ForwardBuckets && fwb, BackwardBuckets && bwb, GroupBounds && gb) :
	MyBaseClass(std::move(fwb), std::move(bwb), std::move(gb))
	{}
	DGCLabeling_Boost(DGCLabeling_Boost &&) = default;
};

class DGCLabeling_Sserialize: public
	DGCLabeling<
	sserialize::Static::Array<
		sserialize::Static::Array<StaticLabelElement>
	>,
	sserialize::Static::Array<
		sserialize::Static::Array<StaticLabelElement>
	>,
	sserialize::Static::Array<int>
	>
{
public:
	using Self = DGCLabeling_Sserialize;
	using MyBaseClass = DGCLabeling<
	sserialize::Static::Array<
		sserialize::Static::Array<StaticLabelElement>
	>,
	sserialize::Static::Array<
		sserialize::Static::Array<StaticLabelElement>
	>,
	sserialize::Static::Array<int>
	>;
	using ForwardBuckets = MyBaseClass::ForwardBuckets;
	using BackwardBuckets = MyBaseClass::BackwardBuckets;
	using GroupBounds = MyBaseClass::GroupBounds;
public:
    static Self fromBinary(sserialize::UByteArrayAdapter data)
    {
        ForwardBuckets fwb;
		BackwardBuckets bwb;
		GroupBounds gb;
		data >> fwb >> bwb >> gb;
		return Self(std::move(fwb), std::move(bwb), std::move(gb));
    }

	
	DGCLabeling_Sserialize(ForwardBuckets && fwb, BackwardBuckets && bwb, GroupBounds && gb) :
	MyBaseClass(std::move(fwb), std::move(bwb), std::move(gb))
	{}
	DGCLabeling_Sserialize(DGCLabeling_Sserialize &&) = default;
};

template<typename TDistanceInfoVector>
class DGCLabeling_Sserialize_Oa: public
	DGCLabeling<
		LabelBuckets,
		LabelBuckets,
		sserialize::Static::Array<int>,
		TDistanceInfoVector
	>
{
public:
	using Self = DGCLabeling_Sserialize_Oa;
	using MyBaseClass = 	DGCLabeling<
		LabelBuckets,
		LabelBuckets,
		sserialize::Static::Array<int>,
		TDistanceInfoVector
	>;
	using ForwardBuckets = typename MyBaseClass::ForwardBuckets;
	using BackwardBuckets = typename MyBaseClass::BackwardBuckets;
	using GroupBounds = typename MyBaseClass::GroupBounds;
public:
    static Self fromBinary(sserialize::UByteArrayAdapter data)
    {
        ForwardBuckets fwb;
		BackwardBuckets bwb;
		GroupBounds gb;
		data >> fwb >> bwb >> gb;
		return Self(std::move(fwb), std::move(bwb), std::move(gb));
    }

	
	DGCLabeling_Sserialize_Oa(ForwardBuckets && fwb, BackwardBuckets && bwb, GroupBounds && gb) :
	MyBaseClass(std::move(fwb), std::move(bwb), std::move(gb))
	{}
	DGCLabeling_Sserialize_Oa(DGCLabeling_Sserialize_Oa &&) = default;
};

template<typename TDistanceInfoVector>
class DGCLabeling_Postgres_Oa: public
	DGCLabeling<
		DBLabels,
		DBLabels,
		std::vector<int>,
		TDistanceInfoVector
	>
{
public:
	using Self = DGCLabeling_Postgres_Oa;
	using MyBaseClass = 	DGCLabeling<
		DBLabels,
		DBLabels,
		std::vector<int>,
		TDistanceInfoVector
	>;
	using ForwardBuckets = typename MyBaseClass::ForwardBuckets;
	using BackwardBuckets = typename MyBaseClass::BackwardBuckets;
	using GroupBounds = typename MyBaseClass::GroupBounds;
public:
	
	DGCLabeling_Postgres_Oa(DataBaseConfig const & cfg, GroupBounds groupBounds = GroupBounds(), std::size_t node_count = 0) :
	MyBaseClass(
				ForwardBuckets(cfg, ForwardBuckets::Type::Forward, node_count),
				BackwardBuckets(cfg, BackwardBuckets::Type::Backward, node_count),
				std::move(groupBounds)
	),
	_con(cfg.connect())
	{
#ifdef GROUP_OF_NODE_FROM_DB
		_con.prepare("get_group", "select g.grp from grp as g where g.pos = $1 limit 1;");
#endif
		if (!groupBounds.size()) {
			pqxx::nontransaction t(_con);
			auto result = t.exec("select bounds.grp, bounds.bound from bounds;");
			this->_groupBounds = std::vector<int>(result.size(), 0);
			for (auto [grp, bound] : result.iter<uint64_t, uint64_t>()) {
				this->_groupBounds.at(grp) = bound;
			}
			this->_numHops = t.query_value<std::size_t>("select stats.group_count from stats limit 1;");
			assert(this->_numHops == this->_groupBounds.size());
		}
	}
	DGCLabeling_Postgres_Oa(Self &&) = default;
protected:
#ifdef GROUP_OF_NODE_FROM_DB
	int _getGroupOfNode(int rank) override {
		pqxx::nontransaction t(_con);
		return t.exec_prepared1("get_group", rank)[0].as<int>();
	}
#endif
private:
	pqxx::connection _con;
};

void help(const char * prog_name) {
    std::cout << "Usage: " << prog_name << "[boost|sv|oav|oam|oaum|oaoah|sqloa|export-v|export-oa|export-sql|export-sql-oa] --query src tgt --dbname <name> -r <num test runs> -i <input binary> -s <seed value> --dca (hop|chhl|ch) --query-file <file> --only-generate-queries --stats -j <num threads> --direct-io --no-mmap --chunked-mmap --advise-random --quiet --checkpoint-size <num>" << std::endl;
}

std::vector<std::pair<int, int>> read_queries(std::filesystem::path const & path) {
	std::ifstream input(path);
	std::vector<std::pair<int, int>> result;
	while(!input.eof() && input.good()) {
		long src, tgt;
		input >> src >> tgt;
		if (input.good()) {
			result.emplace_back(src, tgt);
		}
	}
	return result;
}

void write_queries(std::vector<std::pair<int, int>> const & queries, std::filesystem::path const & path) {
	std::ofstream output(path);
	for(auto [src, tgt] : queries) {
		output << src << '\t' << tgt << '\n';
	}
	output.close();
}

int main(int argc, char *argv[])
{
	std::string op;
	std::string fn;
	std::size_t run_count = 1000;
	std::size_t checkpoint_size = 0;
	DistanceComputationAlgo dca = DistanceComputationAlgo::HopGroups;
	DataBaseConfig db_cfg;
	std::string query_file;
	bool generate_queries = false;
	bool stats = false;
	std::size_t threadCount = 1;
	bool direct_io = false;
	bool no_mmap = false;
	bool chunked_mmap = false;
	bool advise_random = false;
	
	int src = -1;
	int tgt = -1;

	for(int i(1); i < argc; ++i) {
		std::string token(argv[i]);
		if (token == "-i") {
			if (i+1 < argc) {
				fn = argv[i+1];
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "-s") {
			if (i+1 < argc) {
				unsigned int seed = atol(argv[i+1]);
				std::cout << "Using user provided seed value of " << seed << std::endl;
				my_rng.seed(seed);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "-r") {
			if (i+1 < argc) {
				run_count = atol(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "-j") {
			if (i+1 < argc) {
				threadCount = atol(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--checkpoint-size") {
			if (i+1 < argc) {
				checkpoint_size = atol(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--dbhost") {
			if (i+1 < argc) {
				db_cfg.host = std::string(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--dbport") {
			if (i+1 < argc) {
				db_cfg.port = std::string(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--dbname") {
			if (i+1 < argc) {
				db_cfg.name = std::string(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--dbuser") {
			if (i+1 < argc) {
				db_cfg.user = std::string(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--dbpass") {
			if (i+1 < argc) {
				db_cfg.password = std::string(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--query-file") {
			if (i+1 < argc) {
				query_file = std::string(argv[i+1]);
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--query") {
			if (i+2 < argc) {
				src = std::atoi(argv[i+1]);
				tgt = std::atoi(argv[i+2]);
				i += 2;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--dca") {
			if (i+1 < argc) {
				token = std::string(argv[i+1]);
				if (token == "hop" || token == "dgc" || token == "khop") {
					dca = DistanceComputationAlgo::HopGroups;
				}
				else if (token == "chhl") {
					dca = DistanceComputationAlgo::CHHL;
				}
				else if (token == "ch") {
					dca = DistanceComputationAlgo::CH;
				}
				else {
					help(argv[0]);
					std::cerr << "Unkown distance computation algorithm: " << token << std::endl;
					return -1;
				}
				++i;
			}
			else {
				help(argv[0]);
				return -1;
			}
		}
		else if (token == "--only-generate-queries") {
			generate_queries = true;
		}
		else if (token == "--stats") {
			stats = true;
		}
		else if (token == "--direct-io") {
			direct_io = true;
		}
		else if (token == "--no-mmap") {
			no_mmap = true;
		}
		else if (token == "--chunked-mmap") {
			chunked_mmap = false;
		}
		else if (token == "--advise-random") {
			advise_random = true;
		}
		else if (token == "--quiet") {
			quiet = true;
		}
		else if (token == "-h" || token == "--help") {
			help(argv[0]);
			return 0;
		}
		else if (token.size() && token.front() == '-') {
			std::cerr << "Unkown option " << token << std::endl;
			return -1;
		}
		else {
			op = token;
		}
	}
	if (threadCount > 1 && chunked_mmap) {
		std::cerr << "Chunked mmap only supports single threading" << std::endl;
		return -1;
	}
	if (!quiet) {
		std::cout << "Operation: " << op << '\n';
		std::cout << "File: " << fn << '\n';
		std::cout << "dca: ";
		switch (dca) {
			case DistanceComputationAlgo::HopGroups:
				std::cout << "hop";
				break;
			case DistanceComputationAlgo::CHHL:
				std::cout << "chhl";
				break;
			case DistanceComputationAlgo::CH:
				std::cout << "ch";
				break;
			default:
				std::cout << "invalid";
				break;
		}
		std::cout << '\n';
		if (src != -1) {
			std::cout << "src: " << src << '\n';
			std::cout << "tgt: " << tgt << '\n';
		}
	}
	
	auto timedDistance = [&](auto && dgc) {
		auto result = dgc.getTimedDistance(src, tgt);
		std::cout << "dist=" << result.first << '\n';
		std::cout << "time=" << result.second.count() << "us\n";
	};
	auto speedTest = [&](auto && dgc, std::vector<std::pair<int, int>> const & queries) {
		if (!quiet) {
			std::cout << "Running benchmark" << std::endl;
		}
		auto result = dgc.speedTest(queries, dca, threadCount, checkpoint_size);
		std::cout << "mean time per query [us]: " << result.back().second.count()/result.back().first << '\n';
		if (result.size() > 1) {
			std::cout << "count,cumulative time\n";
			for(auto [query_count, comp_time] : result) {
				std::cout << query_count << ',' << comp_time.count() << '\n';
			}
		}
		
		if (!quiet) {
			std::cout << "Queries took on average " << result.back().second.count()/result.back().first << " us." << std::endl;
		}
	};
	auto load_oom_file = [&]() {
		sserialize::UByteArrayAdapter result;
		if (direct_io) {
			result = sserialize::UByteArrayAdapter::open(fn, false, true, 0, 0);
		}
		else if (chunked_mmap) {
			result = sserialize::UByteArrayAdapter::open(fn, false, false, 0);
		}
		else if (no_mmap) {
			result = sserialize::UByteArrayAdapter::open(fn, false, false, 0, 0);
		}
		else {
			result = sserialize::UByteArrayAdapter::open(fn, false, false, std::numeric_limits<std::size_t>::max());
		}
		if (advise_random) {
			result.advice(sserialize::UByteArrayAdapter::AT_RANDOM_READ);
		}
		#if defined(SSERIALIZE_UBA_OPTIONAL_REFCOUNTING)
			if (!quiet) {
				std::cout << "Disabling refcounting of data" << std::endl;
			}
			result.disableRefCounting();
		#endif
		return result;
	};
	auto handle = [&](auto && dgc) {
		if (src != -1) {
			timedDistance(dgc);
		}
		else if (stats) {
			dgc.stats();
		}
		else {
			if (query_file.size()) {
				std::vector<std::pair<int, int>> queries;
				if (std::filesystem::exists(query_file)) {
					queries = read_queries(query_file);
					if (!quiet) {
						std::cout << "Retrieved " << queries.size() << " queries from " << query_file << std::endl;
					}
				}
				else {
					queries = dgc.gen_queries(run_count);
					if (!quiet) {
						std::cout << "Writing " << queries.size() << " queries to " << query_file << std::endl;
					}
					write_queries(queries, query_file);
				}
				if (!generate_queries) {
					speedTest(dgc, queries);
				}
			}
			else {
				speedTest(dgc, dgc.gen_queries(run_count));
			}
		}
	};
	
	if (op == "boost") {
		if (!quiet) {
			std::cout << "Importing data" << std::endl;
		}
		auto dgc = DGCLabeling_Boost::fromBinary(fn);
		if (!quiet) {
			std::cout << "Finished importing data" << std::endl;
		}
		
		handle(dgc);
	}
	else if (op == "sv") {
		if (!quiet) {
			std::cout << "Importing data" << std::endl;
		}
		auto data = load_oom_file();
		auto dgc = DGCLabeling_Sserialize::fromBinary(data);
		if (!quiet) {
			std::cout << "Finished importing data" << std::endl;
		}

		handle(dgc);
	}
	else if (op == "oav") {
		auto data = load_oom_file();
		auto dgc = DGCLabeling_Sserialize_Oa<DistanceInfoVector>::fromBinary(data);
		if (!quiet) {
			std::cout << "Finished importing data" << std::endl;
		}

		handle(dgc);
	}
	else if (op == "oam") {
		auto data = load_oom_file();
		auto dgc = DGCLabeling_Sserialize_Oa<std::map<uint32_t, Distance>>::fromBinary(data);
		if (!quiet) {
			std::cout << "Finished importing data" << std::endl;
		}
		
		handle(dgc);
	}
	else if (op == "oaum") {
		auto data = load_oom_file();
		auto dgc = DGCLabeling_Sserialize_Oa<std::unordered_map<uint32_t, Distance>>::fromBinary(data);
		if (!quiet) {
			std::cout << "Finished importing data" << std::endl;
		}
		
		handle(dgc);
	}
	else if (op == "oaoah") {
		auto data = load_oom_file();
		auto dgc = DGCLabeling_Sserialize_Oa<sserialize::OADHashTable<uint32_t, Distance>>::fromBinary(data);
		if (!quiet) {
			std::cout << "Finished importing data" << std::endl;
		}
		
		handle(dgc);
	}
	else if (op == "sqloa") {
		using DGC = DGCLabeling_Postgres_Oa<sserialize::OADHashTable<uint32_t, Distance>>;
#ifdef SQLOA_GROUP_BOUNDS_FROM_FILE
		auto data = load_oom_file();
		auto dgc_ondisk = DGCLabeling_Sserialize_Oa<sserialize::OADHashTable<uint32_t, Distance>>::fromBinary(data);
		DGC::GroupBounds groupBounds(dgc_ondisk.groubBounds().begin(), dgc_ondisk.groubBounds().end());
		auto dgc = DGC(db_cfg, std::move(groupBounds), dgc_ondisk.numNodes());
#else
		auto dgc = DGC(db_cfg);
#endif
		if (!quiet) {
			std::cout << "Finished importing data" << std::endl;
		}
		
		handle(dgc);
	}
	else if (op == "export-v") {
		std::string outfilename(fn + ".sserialize-vector");
		auto outfile = sserialize::UByteArrayAdapter::createFile(0, outfilename);
		
		std::cout << "Importing data" << std::endl;
		auto dgc = DGCLabeling_Boost::fromBinary(fn);
		std::cout << "Finished importing data" << std::endl;
		std::cout << "Exporting data" << std::endl;
		dgc.sserialize(outfile);
		std::cout << "Finished exporting data" << std::endl;
	}
	else if (op == "export-oa") {
		std::string outfilename(fn + ".sserialize-oa");
		auto outfile = sserialize::UByteArrayAdapter::createFile(0, outfilename);
		
		std::cout << "Importing data" << std::endl;
		auto dgc = DGCLabeling_Boost::fromBinary(fn);
		std::cout << "Finished importing data" << std::endl;
		std::cout << "Exporting data" << std::endl;
		dgc.sserialize_oa(outfile);
		std::cout << "Finished exporting data" << std::endl;
	}
	else if (op == "export-sql" || op == "export-sql-oa") {
		auto my_export = [&](auto && dgc) -> int {
			std::cout << "Exporting data" << std::endl;
			try {
				auto c = db_cfg.connect();
				if (op == "export-sql") {
					dgc.export_postgresql(c);
				}
				else {
					dgc.export_postgresql_oa(c);
				}
			}
			catch(std::exception const & e) {
				std::cout << "Error occured: " << e.what() << std::endl;
				return -1;
			}
			std::cout << "Finished exporting data" << std::endl;
			return 0;
		};
		std::cout << "Importing data" << std::endl;
		if (fn.ends_with("sserialize-oa")) {
			auto data = load_oom_file();
			#if defined(SSERIALIZE_UBA_OPTIONAL_REFCOUNTING)
				std::cout << "Disabling refcounting of data" << std::endl;
				data.disableRefCounting();
			#endif
			auto dgc = DGCLabeling_Sserialize_Oa<sserialize::OADHashTable<uint32_t, Distance>>::fromBinary(data);
			std::cout << "Finished importing data" << std::endl;
			return my_export(dgc);
		}
		else {
			auto dgc = DGCLabeling_Boost::fromBinary(fn);
			std::cout << "Finished importing data" << std::endl;
			return my_export(dgc);
		}
	}
	else {
		std::cerr << "Unkown op " << op << std::endl;
		return -1;
	}
    return 0;
}

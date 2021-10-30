FROM ubuntu:21.04 AS base

ENV DEBIAN_FRONTEND=noninteractive

ARG SRC_DIR="/usr/src"
ENV SRC_DIR="${SRC_DIR}"

ARG IMAGE_DIR="/image"
ENV IMAGE_DIR="${IMAGE_DIR}"

ARG NUM_THREADS=4
ENV NUM_THREADS="${NUM_THREADS:-4}"

ENV POSTGRES_VERSION=13

RUN apt-get update \
    && apt-get -y dist-upgrade \
    && apt-get install wget gnupg2 lsb-core -y \
    && wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add - \
    && echo "deb [ trusted=yes ] https://apt.postgresql.org/pub/repos/apt/ `lsb_release -cs`-pgdg main" | tee /etc/apt/sources.list.d/pgdg.list \
    && apt-get update \
    && apt-get install -y apt-transport-https ca-certificates \
    && apt-get install -y --no-install-recommends \
    postgresql-${POSTGRES_VERSION} \
    postgresql-contrib-${POSTGRES_VERSION} \
    ca-certificates python \
    libcgal-dev libboost-serialization-dev libboost-iostreams-dev

FROM base AS buildenv
RUN apt-get install -y --no-install-recommends \
    postgresql-server-dev-${POSTGRES_VERSION} \
    cmake build-essential ragel protobuf-compiler \
    libicu-dev libcairo-dev libtbb-dev libprotobuf-dev \
    libcrypto++-dev \
    git-core curl

FROM base AS runenv
RUN apt-get install -y --no-install-recommends \
    postgresql-${POSTGRES_VERSION} \
    postgresql-contrib-${POSTGRES_VERSION} \
    libicu67 libcairo2 libtbb2 libprotobuf23 libprotobuf-lite23 \
    libcrypto++  \
    python3-prettytable python3-typeguard python3-psycopg2 \
    sudo cgdb \
    && apt-get clean autoclean \
    && apt-get autoremove --yes \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/

#Install OsmGraphCreator
FROM buildenv AS osmgraphcreator

WORKDIR ${SRC_DIR}
RUN git clone --recursive https://github.com/fmi-alg/OsmGraphCreator OsmGraphCreator
RUN git -C ${SRC_DIR}/OsmGraphCreator fetch \ 
	&& git -C ${SRC_DIR}/OsmGraphCreator checkout 54b8a61035867bc9ef86d8fcfd711ac0b7609a12 \
	&& git -C ${SRC_DIR}/OsmGraphCreator submodule update --init --recursive
WORKDIR ${SRC_DIR}/OsmGraphCreator/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ../
RUN make -j ${NUM_THREADS}

##Copy binaries
RUN mkdir -p ${IMAGE_DIR}/usr/local/bin
RUN cp creator/creator ${IMAGE_DIR}/usr/local/bin/graph-creator 
RUN chmod +x ${IMAGE_DIR}/usr/local/bin/graph-creator
RUN strip ${IMAGE_DIR}/usr/local/bin/graph-creator

##Copy configs
RUN mkdir -p ${IMAGE_DIR}/etc/graph-creator
RUN cp -r ${SRC_DIR}/OsmGraphCreator/data/configs ${IMAGE_DIR}/etc/graph-creator/

#Install CHConstructor
FROM buildenv AS chconstructor

WORKDIR ${SRC_DIR}
RUN git clone --recursive https://github.com/chaot4/ch_constructor.git CHConstructor
RUN git -C ${SRC_DIR}/CHConstructor fetch \ 
	&& git -C ${SRC_DIR}/CHConstructor checkout 05eb8df84b08b2c3f66ec1117d31d44f2e536eca \
	&& git -C ${SRC_DIR}/CHConstructor submodule update --init --recursive
WORKDIR ${SRC_DIR}/CHConstructor/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ../
RUN make -j ${NUM_THREADS}

##Copy binaries
RUN mkdir -p ${IMAGE_DIR}/usr/local/bin
RUN cp ch_constructor ${IMAGE_DIR}/usr/local/bin/ch-constructor
RUN chmod +x ${IMAGE_DIR}/usr/local/bin/ch-constructor
RUN strip ${IMAGE_DIR}/usr/local/bin/ch-constructor

#Install dgc
FROM buildenv AS dgc

WORKDIR ${SRC_DIR}
RUN git clone --recursive -b master https://github.com/fmi-alg/dgc.git
RUN git -C ${SRC_DIR}/dgc fetch \ 
	&& git -C ${SRC_DIR}/dgc checkout 6436dcaf2058f207ce6caa12747ff762526ab076 \
	&& git -C ${SRC_DIR}/dgc submodule update --init --recursive
WORKDIR ${SRC_DIR}/dgc/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DOPTIMIZE_ULTRA=TRUE ../
RUN make -j ${NUM_THREADS}

##Copy binaries
RUN mkdir -p ${IMAGE_DIR}/usr/local/bin
RUN cp dgc-query ${IMAGE_DIR}/usr/local/bin/dgc-query
RUN chmod +x ${IMAGE_DIR}/usr/local/bin/dgc-query
RUN strip ${IMAGE_DIR}/usr/local/bin/dgc-query
RUN cp dgc-create ${IMAGE_DIR}/usr/local/bin/dgc-create
RUN chmod +x ${IMAGE_DIR}/usr/local/bin/dgc-create
RUN strip ${IMAGE_DIR}/usr/local/bin/dgc-create
WORKDIR ${SRC_DIR}/dgc/scripts
# RUN cp bench.py ${IMAGE_DIR}/usr/local/bin/bench.py
# RUN chmod +x ${IMAGE_DIR}/usr/local/bin/bench.py
# RUN cp sql_shopa_func.py ${IMAGE_DIR}/usr/local/bin/sql_shopa_func.py
# RUN chmod +x ${IMAGE_DIR}/usr/local/bin/sql_shopa_func.py

#Final phase
#We simply copy all stuff from ${IMAGE_DIR} of the previous images and add our scripts

FROM runenv AS final

COPY --from=osmgraphcreator ${IMAGE_DIR}/ /
COPY --from=chconstructor ${IMAGE_DIR}/ /
COPY --from=dgc ${IMAGE_DIR}/ /

RUN ldconfig

# Configure PosgtreSQL
RUN mv /etc/postgresql/${POSTGRES_VERSION} /etc/postgresql/current \
  && cd /etc/postgresql/ \
  && ln -s current ${POSTGRES_VERSION}
COPY postgresql.custom.conf.tmpl /etc/postgresql/current/main/
RUN chown -R postgres:postgres /var/lib/postgresql \
 && chown postgres:postgres /etc/postgresql/current/main/postgresql.custom.conf.tmpl \
 && echo "host all all 0.0.0.0/0 md5" >> /etc/postgresql/current/main/pg_hba.conf \
 && echo "host all all ::/0 md5" >> /etc/postgresql/current/main/pg_hba.conf

#Setup users
RUN useradd -m -U dgc

#Setup directories
RUN mkdir "/results" \
    && mkdir "/source"

#Prepare sudo
RUN mkdir -p /etc/sudoers.d

#Copy scripts
COPY drop_caches.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/drop_caches.sh
RUN echo "%dgc     ALL = NOPASSWD: /usr/local/bin/drop_caches.sh" > /etc/sudoers.d/drop_caches

COPY fix_permissions.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/fix_permissions.sh
RUN echo "%dgc     ALL = NOPASSWD: /usr/local/bin/fix_permissions.sh" >> /etc/sudoers.d/fix_perms

COPY sql_shopa_func.py /usr/local/bin/sql_shopa_func.py
RUN chmod +x /usr/local/bin/sql_shopa_func.py

COPY bench.py /usr/local/bin/bench.py
RUN chmod +x /usr/local/bin/bench.py

RUN echo "POSTGRES_VERSION=${POSTGRES_VERSION}" > /run.env.sh 
COPY run.sh /
RUN chmod +x /run.sh

USER root
ENTRYPOINT ["/run.sh"]
CMD []

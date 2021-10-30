#!/bin/bash
service postgresql stop
service postgresql start
echo 1 > /drop_caches

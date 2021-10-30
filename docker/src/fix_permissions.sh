#!/bin/bash
chown -R dgc:dgc /source || exit 1
chown -R dgc:dgc /data || exit 1
chown -R dgc:dgc /results || exit 1
chown -R postgres:postgres /var/lib/postgres || exit 1

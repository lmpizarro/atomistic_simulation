#!/usr/bin/env python

import sys
from pymongo import MongoClient
import time

client = MongoClient()
db = client.pbs_status

sys.path.append("/opt/lmpizarro/python/vtk_examples/SimMatNuc/scripts/libs")

import qsub as qs

# CONFIGURE !!

HOST = "10.128.79.77"
USER = "pizarro"
PASSWORD = "2001llmp"

if __name__ == "__main__":

    host_access = {"user": USER, "host": HOST, "password":PASSWORD}
    remote_commands = [{"in":"qselect -u %s -s R", "out":"", 
        "user":"murina","desc":"R"}, 
                       {"in":"qselect -u %s -s Q", "out":"", 
                           "user":"murina", "desc":"Q"}, {"in":"qstat", "out":"", 
                           "user":"", "desc":"qstat"}]
    outs = qs.qstat_u(host_access, remote_commands)
    qs.qparse_outs (outs)

    for o in outs:
        if o['desc'] == "qstat":
            o['insert_time'] = time.time()
            result = db.qstat.insert_one(o)
            print "id: ", result.inserted_id

    #print (outs)

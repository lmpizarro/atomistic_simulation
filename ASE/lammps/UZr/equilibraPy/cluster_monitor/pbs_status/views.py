from pymongo import MongoClient
import pymongo
import time

client = MongoClient()
db = client.pbs_status

# Create your views here.
from django.shortcuts import render

import sys

def index(request):


    outs = [{"in":"qstat", "out":[{"JobID":12345},{"JobID":12344}], "user":"",
        "desc":"qstat"},
            {"in":"qstat", "out":[{"JobID":12346},{"JobID":12348}], "user":"", "desc":"qstat"}
            ]

    cursor = db.qstat.find().sort('insert_time', pymongo.DESCENDING).limit(1)

    outs = []
    for d in cursor:
        d['hr_time'] = time.ctime(d['insert_time'])
        outs.append(d)


    context = {'qstat_out_list': outs}
 
    print "context ---------------: ", context['qstat_out_list'] 
    return render(request, "index.html", context)

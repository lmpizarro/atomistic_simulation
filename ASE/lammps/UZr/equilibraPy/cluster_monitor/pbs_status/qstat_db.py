from pymongo import MongoClient
import pymongo

client = MongoClient()

db = client.pbs_status


cursor = db.qstat.find().sort('insert_time', pymongo.DESCENDING).limit(1)

outs = []
for d in cursor:
    outs.append(d)



for d in outs:
    for d in d['out']:
        print d['JobID']


print outs    

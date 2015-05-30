#!/usr/bin/env python
import sys
import xmlrpclib
import time

from config import Variables

# if len(sys.argv) < 2 or len(sys.argv) > 3:
#   sys.stderr.write("Usage: NeosClient <xmlfilename | help | queue>\n")
#   sys.exit(1)

neos=xmlrpclib.Server("http://%s:%d" % (Variables.NEOS_HOST, Variables.NEOS_PORT))

# if sys.argv[1] == "help":
#   sys.stdout.write("Help not yet available...\n")

# elif sys.argv[1] == "queue":
#   msg = neos.printQueue()
#   sys.stdout.write(msg)

# else:
# xmlfile = open(sys.argv[1],"r")
xmlfile = open("test.xml", "r")
xml=""
buffer=1
while buffer:
  buffer =  xmlfile.read()
  xml+= buffer
xmlfile.close()

jobNumber = []
password = []
for i in range(10):
  (jN,pd) = neos.submitJob(xml)
  jobNumber.append(jN)
  password.append(pd)
  sys.stdout.write("jobNumber = %d\tpassword = %s\n" % (jobNumber[i],password[i]))

offset=0
status="Waiting"
while status == "Running" or status=="Waiting":
  time.sleep(1)
  # (msg,offset) = neos.getIntermediateResults(jobNumber,password,offset)
  # sys.stdout.write(msg.data)
  status = neos.getJobStatus(jobNumber[9], password[9])

# msg = neos.getFinalResults(jobNumber, password).data
# sys.stdout.write(msg)



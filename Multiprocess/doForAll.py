import subprocess
import multiprocessing
import os
PROC_LIMIT=1
LSTFILE='lst.txt'

def Work(flag,tag):
    subprocess.Popen('echo '+tag,shell=True).wait()


fa=open(LSTFILE)
jobs=[]
i=1

for line in fa:
    tag=line.restrip().split('\t')[0]
    if 1==1:
        p=multiprocessing.Process(target=Work, args=(1,tag))
        p.start()
        jobs.append(p)
        if len(jobs)>=PROC_LIMIT:
            for p in jobs:
                p.join()
            jobs=[]
for p in jobs:
    p.join()
    
    

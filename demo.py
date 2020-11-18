# Generated with SMOP  0.41
from libsmop import *
# demo.m

    a=arange(1,10)
# demo.m:1
    b=concat([2,4,6])
# demo.m:2
    for i in arange(1,length(a)).reshape(-1):
        if sum(a(i) == b) > 0:
            a[i]=0
# demo.m:5
    
    a[a == 0]=[]
# demo.m:8
    for j in arange(1,length(b)).reshape(-1):
        a[a == b(i)]=[]
# demo.m:11
    
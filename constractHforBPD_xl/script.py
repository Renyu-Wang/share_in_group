import os
a=range(1,11)
for n in a:
    command4="./bike2.out %d 3 3 1 1 1 1 1 1"%(n*50)
    os.system(command4)
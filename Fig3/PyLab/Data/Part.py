import numpy as np

class Data(object):
    pass

class Part(object):
    
    def __init__(self):        
        self.data=Data()
    def fields(self,Fields,table='Data'):
        req='SELECT ' + ','.join(Fields) + ' From ' + table  + str(self.sub) + ' WHERE p = ' + str(self.p) + ' ORDER BY t'
        print(req)
        self.fetch(req,Fields)#

    def fetch(self,req,Fields):#
         #print '\r' + req
         self.DB.cursor.execute(req)
         ps=[]
         [ps.append(p) for p in self.DB.cursor]
         for ind,f in zip(range(len(Fields)),Fields):
             #self.data.__setattr__(f,np.array([line[f] for line in self.DB.cursor]))
             #self.data.__setattr__(f,np.array([line[ind] for line in self.DB.cursor]))
             self.data.__setattr__(f,np.array([p[ind] for p in ps]))

    def fieldsgrid(self,Fields,table='Data'):
        req='SELECT ' + ','.join(Fields) + ' From ' + table  + str(self.sub) + ' WHERE nx = ' + str(self.p)
        self.fetch(req,Fields)#




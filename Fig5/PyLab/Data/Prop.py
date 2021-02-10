import numpy as np

class Data(object):
    pass

class Prop(object):
    
    def __init__(self):        
        self.data=Data()
    def fields(self,Fields,table='Prop'):
        req='SELECT ' + ','.join(Fields) + ' From ' + table  + str(self.sub)
        self.fetch(req,Fields)#
        
    def fetch(self,req,Fields):#
         #print '\r' + req
         self.DB.cursor.execute(req)
         ps=[]
         [ps.append(p) for p in self.DB.cursor]
         for ind,f in zip(range(len(Fields)),Fields):
             self.data.__setattr__(f,np.array([p[ind] for p in ps]))
    #def fetch(self,req,Fields):#
    #    self.DB.cursor.execute(req)
    #    for ind,f in zip(range(len(Fields)),Fields):
    #         self.data.__setattr__(f,np.array([line[ind] for line in self.DB.cursor]))
    #for f in Fields:
    #    self.data.__setattr__(f,np.array([line[f] for line in self.DB.cursor]))



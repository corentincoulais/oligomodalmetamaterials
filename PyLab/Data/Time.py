import numpy as np

class Data(object):
    pass

class Time(object):
    def __init__(self):        
        self.data=Data()         
    def fields(self,Fields,table='Data'):
        idx="p"
        if table=='Vert': idx="v"
        req='SELECT ' + ','.join(Fields) + ' From ' + table  + str(self.sub) + ' WHERE t = ' + str(self.t) + ' ORDER BY ' + idx
        print(req)
        self.DB.cursor.execute(req)
        #lines=np.array([line for line in self.DB.cursor])
        ts=[]
        [ts.append(t) for t in self.DB.cursor]
        for ind,f in zip(range(len(Fields)),Fields):
            #self.data.__setattr__(f,np.array([line[ind] for line in self.DB.cursor]))
            self.data.__setattr__(f,np.array([t[ind] for t in ts]))
            #for f in Fields:
            #    self.data.__setattr__(f,np.array([line[f] for line in self.DB.cursor]))
        #self.DB.cursor.close()
    def timelist(self,table='Data'):
        req='SELECT DISTINCT t From ' + table  + str(self.sub)
        #print req
        self.DB.cursor.execute(req)
        Fields=['timelist']
        #for f in Fields:
        #    self.data.__setattr__(f,np.array([line['t'] for line in self.DB.cursor]))
        ts=[]
        [ts.append(t) for t in self.DB.cursor]
        for ind,f in zip(range(len(Fields)),Fields):
            #self.data.__setattr__(f,np.array([line[ind] for line in self.DB.cursor]))
            self.data.__setattr__(f,np.array([t[ind] for t in ts]))
    
    def fetch(self,req,Fields):#
         self.DB.cursor.execute(req)
         ps=[]
         [ps.append(p) for p in self.DB.cursor]
         for ind,f in zip(range(len(Fields)),Fields):
             self.data.__setattr__(f,np.array([p[ind] for p in ps]))

    def selected_fields(self,Fields,table='Data'):
        idx="p"
        if table=='Vert': idx="v"
        req='SELECT ' + ','.join(Fields) + ' From ' + table  + str(self.sub) + ' WHERE t = ' + str(self.t) + ' ORDER BY ' + idx
        print(req)
        self.DB.cursor.execute(req)
        ts=[]
        [ts.append(t) for t in self.DB.cursor]
        for ind,f in zip(range(len(Fields)),Fields):
            self.data.__setattr__(f,np.array([t[ind] for t in ts]))
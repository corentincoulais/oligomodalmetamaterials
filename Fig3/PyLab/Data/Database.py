#from MySQLdb import connect, cursors
import mysql.connector
class DB(object):
    def __init__(self,host='localhost',user='coco',passwd='coco'):	
        self.host=host
        self.user=user
        self.passwd=passwd
    def connect(self,name):
        self.name=name
        self.connection = mysql.connector.connect(user=self.user, password=self.passwd,
                                              host=self.host,database=self.name)
        self.cursor=self.connection.cursor()
#connection = connect(host=self.host,user=self.user,passwd=self.passwd,db=self.name,local_infile = 1)
#self.cursor = connection.cursor(cursors.DictCursor)
#cnx.close()



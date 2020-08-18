#!/usr/bin/python

'''
	Login file should be in the following format (don't forget to make it readable only by you)
	You can find your socket file at /etc/mysql/mysql.conf.d/mysqld.cnf or /etv/httpd/mysql.conf.d/mysqld.cnf
	[client]
	user=yourusername
	password=yourpassowrd
	socket=/var/run/mysqld/mysqld.sock
'''

from MySQLdb import Connect
login_file = "~/.my.cnf"


def connect_to_db():
	conn = Connect(read_default_file=login_file, autocommit=True)
	c = conn.cursor()
	return conn, c

## Install

On Ubuntu 16.04:

    $ sudo apt-get install libpython2.7-dev libpcap-dev cmake python-pip nodejs npm nodejs-legacy libsqlite3-dev sqlite3 pcregrep
    $ sudo -H pip2 install funcparserlib repoze.lru bcrypt flask pysqlite requests

Then:

    $ cd /path/to/susi/web
    $ sudo python2 setup.py install
    $ cd web
    $ . enter.sh
    $ npm install
    $ webpack
    $ sqlite3 data/data.db < schema.sql

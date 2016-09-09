from collections import namedtuple
from contextlib import closing
from functools import wraps
import argparse
import binascii
import datetime
import json
import os
import sqlite3
import time
import traceback
import urllib

import bcrypt
import flask
from flask import \
        Flask, request, session, g, redirect, url_for, abort, render_template, flash
from repoze.lru import lru_cache

import safe_json
import susi

app = Flask(__name__)
app.config.from_object(__name__)
app.config.update(dict(
    # file paths
    DATABASE = os.path.join(app.root_path, 'data/data.db'),
    CAPTURE_DIR = os.path.join(app.root_path, 'data/captures'),
    PCAP_DIR = os.path.join(app.root_path, 'data/pcap'),

    # auth
    USERNAME = 'admin',
    PASSWORD = '$2b$12$OiEBAo5MIA8XOp9riEqn.u7IjQnt0PmVNl9AHCe0Zl3Z.9ByX51t.',
    API_KEY = 'phoorah6ahmiif9Has4A',

    # flask stuff
    SECRET_KEY = 'asdswd1923uasd1u23sajdk',
))
app.config.from_envvar('SUSI_CONFIG', silent=True)

# ================================ #
# ===    app-specific stuff    === #
# ================================ #

def connect_db():
    rv = sqlite3.connect(app.config['DATABASE'])
    rv.row_factory = sqlite3.Row
    return rv

def get_db():
    if not hasattr(g, 'sqlite_db'):
        g.sqlite_db = connect_db()
    return g.sqlite_db

@app.teardown_appcontext
def close_db(error):
    if hasattr(g, 'sqlite_db'):
        g.sqlite_db.close()

def init_db():
    with closing(connect_db()) as db:
        with app.open_resource('schema.sql', mode='r') as f:
            db.cursor().executescript(f.read())
        db.commit()

def check_login(user, pw):
    pwhash = app.config['PASSWORD']
    return user == app.config['USERNAME'] and bcrypt.hashpw(pw.encode('utf-8'), pwhash) == pwhash

def correct_api_key():
    return request.headers.get('x-api-key') == app.config['API_KEY']

def get_user():
    if 'user' in session:
        return session.get('user')
    if correct_api_key():
        return 'api'
    return None

def generate_token(n):
    return binascii.hexlify(os.urandom(n))

def fetch_all(sql, args=[]):
    cur = get_db().execute(sql, args)
    return [dict(**row) for row in cur.fetchall()]

def fetch_single(sql, args):
    cur = get_db().execute(sql, args)
    res = cur.fetchall()
    assert len(res) == 1, 'selected object does not exist'
    return dict(**res[0])

def delete_single(sql, args):
    db = get_db()
    cur = db.execute(sql, args)
    db.commit()

def insert_and_fetch(table, sql, args):
    db = get_db()
    cur = db.execute(sql, args)
    db.commit()
    return fetch_single('select * from %s where id = ?' % table, [cur.lastrowid])

@lru_cache(64)
def load_collection(dir_):
    return susi.Collection(dir_)

def format_date(ts):
    return datetime.datetime.fromtimestamp(ts/1000000.0).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]

def paginate(result_list, pagination_opts):
    from_ = pagination_opts['from']
    to = pagination_opts['to']
    assert is_int(from_)
    assert to is None or is_int(to)
    if to is None:
        return result_list[from_:], len(result_list)
    else:
        return result_list[from_:to], len(result_list)

def get_collection_from_id(collection_id):
    collection = fetch_single('select * from collections where id=?', [collection_id])
    #collection['first_packet_time'] = format_date(collection['first_packet_time'])
    return collection

def is_int(x):
    return isinstance(x, (int, long))

# ================================ #
# ===        middleware        === #
# ================================ #

def no_login_required(f):
    f.no_login = True
    return f

Success = namedtuple('Success', ('value',))
Failure = namedtuple('Failure', ('value',))

def format_error(e):
    return type(e).__name__ + ': ' + str(e)

def json_route(f):
    @wraps(f)
    def json_wrapper(*args, **kwargs):
        try:
            if request.data:
                request.data = json.loads(request.data)
            if request.query_string:
                request.data = json.loads(urllib.unquote(request.query_string))
        except ValueError:
            traceback.print_exc()
            return flask.Response('invalid json data', status=400)
        try:
            res = f(*args, **kwargs)
        except Exception, e:
            traceback.print_exc()
            res = Failure(format_error(e) if app.config['DEBUG'] else None)
        res = {
            'success': isinstance(res, Success),
            ('result' if isinstance(res, Success) else 'error'): res.value
        }
        encoded = safe_json.dumps(res)
        return flask.Response(encoded, headers={
            'Content-Type': 'application/json'
        })
    return json_wrapper

@app.before_request
def no_login():
    login_required = True
    if request.endpoint:
        login_required = (
            not getattr(app.view_functions[request.endpoint], 'no_login', False)
            and request.endpoint != 'static'
        )
    if login_required and not get_user():
        return redirect(url_for('login'))

@app.before_request
def xsrf_protection():
    token = session.get('xsrf_token')
    if token is None:
        session['xsrf_token'] = token = generate_token(16)

    if request.method == 'GET':
        return
    if correct_api_key():
        return

    user_token = request.form.get('xsrf_token')
    if request.headers.get('x-xsrf-token'):
        user_token = request.headers.get('x-xsrf-token')

    if token != user_token:
        return 'invalid xsrf token'

# ================================ #
# ===          routes          === #
# ================================ #

@app.route('/login', methods=['GET', 'POST'])
@no_login_required
def login():
    error=None
    if request.method == 'POST':
        if check_login(request.form['username'], request.form['password']):
            session['user'] = 'admin'
            return redirect(url_for('index'))
        else:
            error='Invalid username or password'
    return render_template('login.html', error=error)

@app.route('/logout', methods=(['GET', 'POST'] if app.config['DEBUG'] else ['POST']))
@no_login_required
@json_route
def logout():
    session.pop('user')
    return Success(None)

@app.route('/')
def index():
    return render_template('index.html')

def insert_collection(obj):
    assert not '..' in obj['directory']
    directory = obj['directory']

    sql = 'insert into collections (name, directory) values (?, ?)'
    args = [obj['name'], directory]
    return insert_and_fetch('captures', sql, args)

@app.route('/api/collections', methods=['GET', 'POST'])
@json_route
def api_collections():
    if request.method == 'GET':
        collections = fetch_all('select id from collections order by created_at desc')
        return Success([get_collection_from_id(coll['id']) for coll in collections])
    elif request.method == 'POST':
        return Success(insert_collection(request.data))

@app.route('/api/documents', methods=['GET'])
@json_route
def api_get_documents():
    ids = request.data['collection_ids']
    filter_ = request.data['filter']
    pagination = request.data['result_range']
    #sort_by = request.data['sort_by']
    assert isinstance(ids, list) and all(map(is_int, ids))
    assert isinstance(filter_, unicode)

    print "query", filter_
    all_results = []
    query_time = 0
    for capture_id in ids:
        coll = get_collection_from_id(capture_id)
        index = load_collection(coll['directory'])
        time_start = time.time()
        try:
            query_results = index.lookup(filter_, k=10)
        except Exception, e:
            if isinstance(e, (susi.LexerError, susi.NoParseError)):
                return Failure('Error parsing query: ' + str(e))
            else:
                raise
        query_time += time.time() - time_start
        all_results += [{'id': docid,
                         'weight': weight,
                         'snippet': snippet,
                         'collection': coll}
                        for docid, weight, snippet in query_results]

    all_results = sorted(all_results, key=lambda doc: (-int(doc['weight']), doc['id']))
    query_results, result_count = paginate(all_results, pagination)

    return Success({
        'query_time': query_time,
        'matches': query_results,
        'result_count': result_count,
    })


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--port', default=5000, type=int, help='Port')
    p.add_argument('--host', default='127.0.0.1', type=str, help='Hostname')
    p.add_argument('--debug', action='store_true', help='Debug mode')
    args = p.parse_args()

    app.run(debug=args.debug, host=args.host, port=args.port)

var MINI = require('minified');
var $ = MINI.$, _ = MINI._;
var xsrf_token = $('meta[name="xsrf_token"]').get('content');

function request(method, path, data) {
  if (typeof data !== 'undefined')
    data = JSON.stringify(data);
  return $.request(method, path, data, {
    'headers': {
      'X-XSRF-Token': xsrf_token
    }
  }).then(function (data) {
    return JSON.parse(data);
  });
}
function get(path, data) { return request('GET', path, data); }
function post(path, data) { return request('POST', path, data); }
function put(path, data) { return request('PUT', path, data); }
function delete_(path, data) { return request('DELETE', path, data); }

var Spinner = React.createClass({
  render: function() {
    return <span className="glyphicon glyphicon-refresh glyphicon-refresh-animate"></span>;
  },
});

var StreamInfoRow = React.createClass({
  render: function() {
    return (
      <div className={"row info" + (this.props.toggled ? ' toggled' : '')} onClick={this.props.onClick}>
        <div className="col-xs-4">{this.props.stream.capture.name}</div>
        <div className="col-xs-3">{this.props.stream.capture.host}:{this.props.stream.capture.port}</div>
        <div className="col-xs-4">{this.props.stream.timestamp}</div>
        <div className="col-xs-1">{this.props.stream.id}</div>
      </div>
    );
  }
});

var AsciiStreamData = React.createClass({
  render: function() {
    var rows = [];
    var placeholder = <div className="col-xs-1"></div>;
    this.props.stream.chunks.forEach(function(chunk, index) {
      var client_to_server = (index % 2 == 0);
      if (chunk.length === 0)
        return;
      rows.push(
        <div key={index} className={'row ' + (client_to_server ? 'ascii-client' : 'ascii-server')}>
          {client_to_server ? '' : placeholder}
          <div className="col-xs-11 ascii-data">
            <pre>
              {chunk}
            </pre>
          </div>
          {client_to_server ? placeholder : ''}
        </div>
      );
    });
    return (
      <div className="container ascii-viewer">
        {rows}
      </div>
    );
  }
});

// TODO implement hex viewer
var StreamDataViewer = React.createClass({
  getInitialState: function() {
    return {
      stream_details: null,
    };
  },

  componentDidMount: function() {
    get('/api/stream/' + this.props.stream.capture.id + '/' + this.props.stream.id)
    .then(function(result) {
      console.log(result);
      if (result.success) {
        this.setState({ stream_details: result.result });
      }
    }.bind(this));
  },

  render: function() {
    return (
      <div className="row">
        <div className="col-xs-12">
          {this.state.stream_details !== null
            ? <AsciiStreamData stream={this.state.stream_details} />
            : <Spinner />}
        </div>
      </div>
    );
  }
});

var StreamList = React.createClass({
  getInitialState: function() {
    return {
      shown_streams: {},
    };
  },

  getStreamKey: function(stream) {
    return '' + stream.capture.id + '/' + stream.id;
  },

  onRowClick: function(stream) {
    var key = this.getStreamKey(stream);
    var shown_streams = _.extend({}, this.state.shown_streams);
    if (shown_streams[key])
      delete shown_streams[key];
    else
      shown_streams[key] = true;
    this.setState({
      shown_streams: shown_streams,
    });
  },

  render: function() {
    var rows = [];
    this.props.streams.forEach(function(stream) {
      var key = this.getStreamKey(stream);
      rows.push(
        <StreamInfoRow
          key={key} stream={stream} onClick={this.onRowClick.bind(this, stream)}
          toggled={this.state.shown_streams[key] === true}
          />);
      if (this.state.shown_streams[key]) {
        rows.push(<StreamDataViewer key={key + '_dropdown'} stream={stream} />);
      }
    }.bind(this));
    return (
      <div className="container stream-list">
        <div className="row">
          <div className="col-xs-4 th">Capture</div>
          <div className="col-xs-3 th">Host:Port</div>
          <div className="col-xs-4 th">Time</div>
          <div className="col-xs-1 th">SID</div>
        </div>
        {rows}
      </div>
    );
  }
});

var FilterableStreamList = React.createClass({
  getInitialState: function() {
    return {
      streams: [],
      error: null,
      filter: this.props.filter,
      loading: true,
      query_time: null,
    };
  },

  load: function() {
    this.setState({ loading: true });

    get('/api/streams', {
      capture_ids: this.props.capture_ids,
      filter: this.state.filter,
      result_range: { from: 0, to: null },
    }).then(function(result) {
      if (result.success) {
        this.setState({
          streams: result.result.matches,
          query_time: result.result.query_time,
          error: null,
        });
      } else {
        this.setState({
          streams: this.state.streams,
          query_time: null,
          error: result.error,
        });
      }
      this.setState({ loading: false });
    }.bind(this));
  },

  componentDidMount: function() {
    this.load();
  },

  handleFilterSubmit: function(evt) {
    evt.preventDefault();
    this.setState({ filter: this.refs.filter.value }, function() {
      this.load();
    }.bind(this));
  },

  render: function() {
    var t = this.state.query_time;
    var t = t !== null ? Math.round(t*100000)/100 : null;
    return (
      <div>
        <form onSubmit={this.handleFilterSubmit}>
          <div className="input-group">
            <input
              className="form-control"
              id="filter" ref="filter" type="text" defaultValue={this.props.filter}
              autocomplete="new-password"
              />
            <span className="input-group-btn">
              <input type="submit" className="btn btn-default" value="Filter" />
            </span>
          </div>

          <div className="bigger">
            {this.state.loading
              ? <Spinner />
              : ''}
            {this.state.error !== null
              ? <span className="text-danger">{this.state.error}</span>
              : <span>{this.state.streams.length} streams match</span> }
            {t !== null
              ? <span> (query took {t} ms)</span>
              : ''}
          </div>
        </form>
        <StreamList streams={this.state.streams} />
      </div>
    );
  }
});

var initial_filter = '"HTTP" and not "requests"'
ReactDOM.render(
  <FilterableStreamList capture_ids={[19]} filter={initial_filter} />,
  document.getElementById('index-page')
);

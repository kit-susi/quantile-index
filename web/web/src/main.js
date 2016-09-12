import React from 'react';
import ReactDOM from 'react-dom';
import ReactPaginate from 'react-paginate';
import SelectBox from 'react-select-box';
import { $, _ } from 'minified';

let xsrf_token = $('meta[name="xsrf_token"]').get('content');

function request(method, path, data) {
  if (typeof data !== 'undefined')
    data = JSON.stringify(data);
  return $.request(method, path, data, {
    'headers': { 'X-XSRF-Token': xsrf_token }
  }).then((data) => JSON.parse(data));
}

function get(path, data) { return request('GET', path, data); }
function post(path, data) { return request('POST', path, data); }
function put(path, data) { return request('PUT', path, data); }
function delete_(path, data) { return request('DELETE', path, data); }

function arraysEqual(a, b) {
  if (a.length != b.length)
    return false;
  a.every((x, i) => x === b[i]);
}

class Spinner extends React.Component {
  render() {
    return (
      <span
        className="glyphicon glyphicon-refresh glyphicon-refresh-animate"
        style={this.props.style}>
      </span>
    );
  }
}

class DocInfoRow extends React.Component {
  render() {
    let doc = this.props.doc;
    return (
      <div className={"row info" + (this.props.toggled ? ' toggled' : '')} onClick={this.props.onClick}>
        <div className="col-xs-3">{doc.collection.name}</div>
        <div className="col-xs-2">{doc.id}</div>
        <div className="col-xs-3">{doc.snippet.split('\n')[0]}</div>
        <div className="col-xs-2">{doc.weight}</div>
        <div className="col-xs-2">{doc.snippet.length}</div>
      </div>
    );
  }
}

class DocDataViewer extends React.Component {
  constructor(props, context) {
    super(props, context);
  }

  render() {
    return (
      <div className="row viewer">
        <div className="col-xs-12">
          <pre>{this.props.doc.snippet}</pre>
        </div>
      </div>
    );
  }
}

class DocList extends React.Component {
  constructor(props, context) {
    super(props, context);
    this.state = {
      toggled_docs: {}
    };
  }

  getDocKey(doc) {
    return '' + doc.collection.id + '/' + doc.id;
  }

  onRowClick(doc) {
    let key = this.getDocKey(doc);
    let toggled_docs = _.extend({}, this.state.toggled_docs);
    if (toggled_docs[key])
      delete toggled_docs[key];
    else
      toggled_docs[key] = true;
    this.setState({
      toggled_docs: toggled_docs,
    });
  }

  componentWillReceiveProps(props) {
    let new_doc_ids = props.docs.map((doc) => doc.id);
    let old_doc_ids = this.props.docs.map((doc) => doc.id);
    if (!arraysEqual(new_doc_ids, old_doc_ids))
      this.setState({ toggled_docs: {}});
  }

  render() {
    let rows = [];
    this.props.docs.forEach((doc) => {
      let key = this.getDocKey(doc);
      rows.push(
        <DocInfoRow
          key={key} doc={doc} onClick={() => this.onRowClick(doc)}
          toggled={this.state.toggled_docs[key] === true}
          />);
      if (this.state.toggled_docs[key]) {
        rows.push(<DocDataViewer key={key + '_dropdown'} doc={doc} />);
      }
    });
    let columns = [
      [3, 'collection_name', 'Collection'],
      [2, 'doc_id', 'Doc ID'],
      [3, 'title', 'Title'],
      [2, 'weight', 'TF'],
      [2, 'size', 'Size'],
    ];

    return (
      <div className="container doc-list">
        <div className="row">
          {columns.map(
            ([width, column_id, column_text]) => (
              <div key={'col_' + column_id} className={'col-xs-' + width + ' th'}>{column_text}</div>
            ))}
        </div>
        {rows}
      </div>
    );
  }
}

class FilterableDocList extends React.Component {
  constructor(props, context) {
    super(props, context);
    this.state = {
      docs: [],
      error: null,
      filter: this.props.initial_filter,
      loading: true,
      query_time: null,
      offset: 0,
      collection_ids: this.props.collection_ids,
    };
  }

  load() {
    this.setState({ loading: true });

    get('/api/documents', {
      collection_ids: this.state.collection_ids,
      filter: this.state.filter,
      result_range: {
        from: this.state.offset,
        to: this.state.offset + this.props.docs_per_page
      },
    }).then((result) => {
      if (result.success) {
        this.setState({
          docs: result.result.matches,
          query_time: result.result.query_time,
          error: null,
          result_count: result.result.result_count,
        });
      } else {
        this.setState({
          docs: this.state.docs,
          query_time: null,
          error: result.error,
          result_count: 0,
        });
      }
      this.setState({ loading: false });
    });
  }

  componentWillReceiveProps(props) {
    if (!arraysEqual(this.state.collection_ids, props.collection_ids))
      this.setState({ collection_ids: props.collection_ids }, () => this.load())
  }

  componentDidMount() {
    this.load();
    this.refs.filter.focus();
  }

  handleFilterSubmit(evt) {
    evt.preventDefault();
    this.setState({
      filter: this.refs.filter.value,
      offset: 0
    }, () => this.load());
  }

  handlePageClick(data) {
    let selected = data.selected;
    let offset = Math.ceil(selected * this.props.docs_per_page);

    this.setState({offset: offset}, () => this.load());
  };

  render() {
    let t = this.state.query_time;
    t = t !== null ? Math.round(t*100000)/100 : null;
    let to = -1 + Math.min(
        this.state.offset + this.props.docs_per_page,
        this.state.result_count);
    return (
      <div>
        <form onSubmit={(evt) => this.handleFilterSubmit(evt)}>
          <div className="input-group">
            <input
              className="form-control"
              id="filter" ref="filter" type="text" defaultValue={this.state.filter}
              autoComplete="new-password"
              />
            <span className="input-group-btn">
              <input type="submit" className="btn btn-default" value="Filter" />
            </span>
          </div>

          <div className="bigger">
            {this.state.error !== null
              ? <span className="text-danger">{this.state.error}</span>
              : <span>Showing results {this.state.offset} to {to} of <b>
                  {this.state.result_count}</b> matches</span>}
            {t !== null
              ? <span> (query took <b>{t} ms</b>)</span>
              : ''}
            {this.state.loading
              ? <Spinner style={{marginLeft: '10px'}} />
              : ''}
          </div>
        </form>
        <DocList
          docs={this.state.docs}
          />
        { this.state.result_count > this.props.docs_per_page
          ?  <ReactPaginate
                previousLabel={"previous"}
                nextLabel={"next"}
                breakLabel={<a href="">...</a>}
                pageNum={Math.ceil(this.state.result_count / this.props.docs_per_page)}
                marginPagesDisplayed={2}
                pageRangeDisplayed={5}
                clickCallback={(data) => this.handlePageClick(data)}
                containerClassName={"pagination"}
                subContainerClassName={"pages pagination"}
                activeClassName={"active"} />
          : ''}
      </div>
    );
  }
}

class CollectionSelectBox extends React.Component {
  constructor(props, context) {
    super(props, context);
    this.state = {
      collections: [],
      selected_collections: [],
      loading: true,
    };
  }

  collectionListChanged(selected_collections) {
    this.setState({
      selected_collections: selected_collections
    });
    this.props.onCollectionListChange(selected_collections);
  }

  componentDidMount() {
    get('/api/collections')
    .then((result) => {
      if (result.success) {
        this.setState({
          loading: false,
          collections: result.result,
        });
      }
    });
  }

  render() {
    if (this.state.loading)
      return <Spinner />;

    let options = this.state.collections.map((collection) =>
      (
        <option key={collection['id']} value={collection['id']}>
          {collection['name']}
        </option>
      ));
    return (
      <SelectBox
          label="Collections"
          value={this.state.selected_collections}
          multiple={true}
          onChange={(selection) => this.collectionListChanged(selection)}
          >
        {options}
      </SelectBox>
    );
  }
}

class Susi extends React.Component {
  constructor(props, context) {
    super(props, context);
    this.state = {
      collection_ids: [],
      docs_per_page: 20,
      initial_filter: '',
      loading: true,
    };
  }

  collectionListChanged(collection_ids) {
    this.setState({ collection_ids: collection_ids });
  }

  render() {
    return (
      <div>
        <div style={{marginBottom: '10px'}}>
          <b>Select Collections:</b>{' '}
          <CollectionSelectBox
            onCollectionListChange={(collection_ids) => this.collectionListChanged(collection_ids) }/>
        </div>
        <FilterableDocList
          collection_ids={this.state.collection_ids}
          initial_filter={this.state.initial_filter}
          docs_per_page={20}
          />
      </div>
    );
  }
}

ReactDOM.render(
  <Susi />,
  document.getElementById('index-page')
);

drop table if exists collections;
create table collections (
  id integer primary key autoincrement,
  name text,
  directory text,
  created_at timestamp default current_timestamp
);

#!/bin/bash
#
# ./sdsl.sh pull [branch]
# Pull from the given branch in git@github.com:niklasb/susi-sdsl-lite.git
# into subtree /external/sdsl-lite/
# (branch defaults to master)
#
# ./sdsl.sh push [branch]
# Push changes in subtree /external/sdsl-lite/ to
# the given branch in git@github.com:niklasb/susi-sdsl-lite.git
# (branch defaults to master)

set -e
if [[ ( "$1" != "pull" && "$1" != "push" ) || $# > 2 ]]; then
  echo >&2 "Usage: ./sdsl.sh (pull | push) [branch]"
  exit 1
fi

action="$1"
upstream_branch="$2"

if [ ! -d external/sdsl-lite ]; then
  echo >&2 "external/sdsl-lite does not exist"
  exit 1
fi

cur_branch="`git rev-parse --abbrev-ref HEAD`"
if [[ "$cur_branch" != "master" ]]; then
  echo >&2 -n "You are not in master. Are you sure you still want to do this?"
  read -p " [y/N] " -n 1 -r
  echo
  [[ $REPLY =~ ^[Yy]$ ]] || exit 0
fi

# Add sdsl remote if it does not exist
if ! git remote | egrep '^sdsl$' &>/dev/null; then
  set -x
  git remote add sdsl git@github.com:niklasb/susi-sdsl-lite.git
  set +x
fi

set -x
git fetch sdsl
set +x

if [[ "$action" == "pull" ]]; then
  set -x
  git merge -s subtree --squash --no-commit --allow-unrelated-histories \
    sdsl/master
  git commit -m 'Merge upstream sdsl-lite'
  set +x
elif [[ "$action" == "push" ]]; then
  # Mirror sdsl/master into tmp-sdsl-master branch, merge subtree into
  # it, then push it to upstream master
  set -x
  git checkout -b tmp-sdsl-master "sdsl/$upstream_branch"
  git merge -s subtree --squash --no-commit --allow-unrelated-histories \
    "$cur_branch"
  git commit -m 'Merge sdsl-lite subtree from susi into upstream'
  git push sdsl tmp-sdsl-master:"$upstream_branch"
  git checkout "$cur_branch"
  git branch -D tmp-sdsl-master
  set +x
fi

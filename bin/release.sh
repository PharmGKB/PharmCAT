#!/bin/bash
#
# Performs a release on PharmGKB/PharmCAT.
#
#
set -e
set -u
set -o pipefail


branch=$(git rev-parse --abbrev-ref HEAD)
if [[ $branch != 'development' ]]
then
  echo "Must run on development branch."
  exit 1
fi

diffs=$(git diff | wc -l)
if [[ $diffs -gt 0 ]]
then
  echo "There are uncommitted changes."
  exit 1
fi

commits=$(git rev-list HEAD...origin/development | wc -l)
if [[ $commits -gt 0 ]]
then
  echo "This branch is ${commits} behind origin/development."
  exit 1
fi

# update disclaimer.hbs
script_dir="$( dirname -- "$BASH_SOURCE"; )";
"${script_dir}/update_disclaimer.js"

git_config=$(git config --list)
safe_crlf=""
if [[ "${git_config}" == *"core.safecrlf"* ]]
then
  safe_crlf=$(git config core.safecrlf)
  if [[ $safe_crlf == 'warn' ]]
  then
    git config set core.safecrlf false
  fi
fi

diffs=$(git diff | grep -c "${script_dir}/disclaimer.hbs")
if [[ $diffs -gt 0 ]]
then
  git add src/main/resources/org/pharmgkb/pharmcat/reporter/disclaimer.hbs
  git commit -m "chore: update disclaimer"
  git push
fi

if [[ $safe_crlf == 'warn' ]]
then
  git config set core.safecrlf warn
fi


# do release
echo "Creating release..."
yarn release


# update repo
echo ""
echo "Updating repo..."
git fetch


# update main branch
echo ""
echo "Updating main branch..."
git checkout main
git pull
git rebase development
git push

# switching back to development
git checkout development


echo ""
echo "Done."

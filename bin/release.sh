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

# update disclaimers.hbs
echo "Updating disclaimer"
script_dir="$( dirname -- "$BASH_SOURCE"; )";
"${script_dir}/update_disclaimer.js"

# this reverts files with only EOL changes
git stash
git stash pop

diffs=$(git diff --name-only | wc -l)
if [[ $diffs -gt 0 ]]
then
  echo "Committing disclaimer"
  git add src/main/resources/org/pharmgkb/pharmcat/reporter/disclaimers.hbs
  git commit -m "chore: update disclaimer"
  git push
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

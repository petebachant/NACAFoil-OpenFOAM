#!/usr/bin/env bash

cp figures/*.png app/figures
cp processed/all-simulated.csv app/processed
cp processed/NACA0012_6e6_Ladson_180grit.csv app/processed
cp requirements.txt app
cp notebook.py app/app.py
cp layouts/notebook.grid.json app/layouts

# Commit and push if anything has changed
cd app && git add .
DIFF=$(git diff --staged)
if [ "$DIFF" ]
then
  echo "Changes detected, committing and pushing"
  git commit -m "Update app" && git push
else
  echo "No changes detected, skipping commit and push"
fi

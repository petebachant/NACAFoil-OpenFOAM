#!/usr/bin/env bash

cp figures/*.png app/figures
cp processed/all-simulated.csv app/processed
cp processed/NACA0012_6e6_Ladson_180grit.csv app/processed
cp requirements.txt app
cp notebook.py app/app.py
cp layouts/notebook.grid.json app/layouts

cd app && git add . && git commit -m "Update app" && git push

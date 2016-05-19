@echo off
python "fit_from_database.py"
python "aggregate_summary.py"
python "fit_to_summary.py"
python "build_summary_activation.py"
pause
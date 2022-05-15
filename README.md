## Pande-mic DASHBOARD

this is repository contains a Flask application dashboard designed to show the results of the analysis done by the monitoring pipeline [available here](https://github.com/fmfi-compbio/pande-mic)


To run the dashboard app, activate the conda environment you have created for the monitoring pipeline and run
```
python3 ./dashboard.py \
 --summary_path /<path to the summary dir created by the pipeline>/summary/ \
 --config_path /<path to the pipeline config dir>/config/ \
 --run_name "<anything>"
```

Options for the script are:
```
usage: dashboard.py [-h] --summary_path SUMMARY_PATH --config_path CONFIG_PATH --run_name RUN_NAME
                    [--smoothing_window_size SMOOTHING_WINDOW_SIZE]
                    [--smoothing_window_step SMOOTHING_WINDOW_STEP] [--low_cutoff LOW_CUTOFF]
                    [--high_cutoff HIGH_CUTOFF] [--coverage_cutoff] [--port PORT] [--debug]

Flask-app dashboard for pande-mic.

optional arguments:
  -h, --help            show this help message and exit
  --summary_path SUMMARY_PATH
                        full path to summary directory
  --config_path CONFIG_PATH
                        full path to config directory
  --run_name RUN_NAME   run name or ID
  --smoothing_window_size SMOOTHING_WINDOW_SIZE
                        window size used for smoothing of coverage graphs
  --smoothing_window_step SMOOTHING_WINDOW_STEP
                        step between smoothing windows for coverage graphs (calculating madians for each
                        position might be slow)
  --low_cutoff LOW_CUTOFF
                        coverage value that is considered too low
  --high_cutoff HIGH_CUTOFF
                        coverage value that is considered too high
  --coverage_cutoff     cut values higher than high_cutoff from the graph
  --port PORT
  --debug               run app in debug mode

```

To test the dashboard app on the sample output from the pipeline, run

```
python3 ./dashboard.py --summary_path ./sample_pipeline_output/SARS-CoV-2/ --config_path ./sample_pipeline_output/SARS-CoV-2-config/ --run_name "SARS-CoV-2"

```
in case of any errors, you can try:

```
python3 ./dashboard.py --summary_path ./sample_pipeline_output/SARS-CoV-2/ --config_path ./sample_pipeline_output/SARS-CoV-2-config/ --run_name "SARS-CoV-2" --debug
```
this will run the dashboard app with the Flask development server
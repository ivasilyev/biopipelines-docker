#!/usr/bin/env bash

# Create the Jupyter Lab server instance
jupyter lab --ip=0.0.0.0 --port=31522 --no-browser --NotebookApp.token=aJbR9v

# Automatically restart the server
bash $(realpath "$0")

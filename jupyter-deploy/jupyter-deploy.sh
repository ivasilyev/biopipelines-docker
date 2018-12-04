#!/usr/bin/env bash

PORT=31522
TOKEN="aJbR9v"

# Create the Jupyter Lab server instance
jupyter lab --ip=0.0.0.0 --port=${PORT} --no-browser --NotebookApp.token=${TOKEN}

# Automatically restart the server
bash $(realpath "$0")

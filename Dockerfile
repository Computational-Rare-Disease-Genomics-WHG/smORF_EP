# Use the python 3.11.4-slim-bullseye image
FROM python:3.11.4-slim-bullseye

# Install required system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Add python build tools
RUN python3 -m pip install --upgrade pip setuptools wheel

# Set the working directory to /smorf-ep
WORKDIR /smorf-ep

# Copy the current directory contents into the container at /smorf-ep
COPY . /smorf-ep

# Run the command to install dependencies
RUN python3 setup.py install

# The following is a Dockerfile 
# to run smORF-EP on a Docker container

# Use the python 3.11 image slim-buster
FROM python:3.11-slim-buster

# Install build dependencies
RUN apt-get update && apt-get upgrade

# Add python build tools
RUN python3 -m pip install --upgrade pip setuptools wheel

# Set the working directory to /app
WORKDIR /smorf-ep

# Copy the current directory contents into the container at /app
ADD . /smorf-ep

# Run the command to install dependencies
RUN python3 setup.py install

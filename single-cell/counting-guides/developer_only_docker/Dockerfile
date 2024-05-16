# Use the official Ubuntu base image
FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Replace default repositories with a different mirror
RUN sed -i 's|http://ports.ubuntu.com/ubuntu-ports|http://archive.ubuntu.com/ubuntu|g' /etc/apt/sources.list

# Install basic dependencies to update the repository keys
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    gnupg2 \
    && apt-get clean

# Add the keyring for Ubuntu repositories
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 3B4FE6ACC0B21F32 && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 40976EAF437D05B5 && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 46181433FBB75451

# Update and install the rest of the dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    bzip2 \
    curl \
    git \
    && apt-get clean

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
    chmod +x /miniconda.sh && \
    /miniconda.sh -b -p /opt/conda && \
    rm /miniconda.sh

# Add conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Update conda
RUN conda update -n base conda -y

# Install conda-libmamba-solver
RUN conda install -n base conda-libmamba-solver -y

# Set solver to libmamba
RUN conda config --set solver libmamba

# Install Nextflow
RUN conda install -c bioconda nextflow -y

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]

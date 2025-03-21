FROM condaforge/mambaforge:latest

# Copy environment file
COPY vsgseq2.yml /tmp/vsgseq2.yml

# Install dependencies using mamba
RUN mamba env create -f /tmp/vsgseq2.yml && \
    mamba clean -afy

# Add conda environment to PATH
ENV PATH /opt/conda/envs/vsgseq2-env/bin:$PATH

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"]

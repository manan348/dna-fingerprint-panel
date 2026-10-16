# DNA Fingerprint Panel - Reproducible Container
# Build: docker build -t fingerprint-panel:1.0.0 .
# Run: docker run -it --rm -v $(pwd):/work fingerprint-panel:1.0.0

FROM mambaorg/micromamba:1.5.0-bookworm-slim

LABEL maintainer="Abdul Manan <abdul.manan@neurocrop.com>"
LABEL description="DNA Fingerprint Panel for Varietal Identification"
LABEL version="1.0.0"

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    curl \
    ca-certificates \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

USER mamba

# Copy environment file and create conda environment
COPY --chown=mamba:mamba environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

# Set environment
ENV PATH="/opt/conda/bin:${PATH}"
ENV PYTHONPATH="/work/src:${PYTHONPATH}"

# Copy project files
WORKDIR /work
COPY --chown=mamba:mamba . /work

# Install the package in development mode
RUN pip install -e .

# Default command
CMD ["fingerprint", "--help"]
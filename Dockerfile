FROM mambaorg/micromamba:2-alpine3.22

ARG ODT_REF=main
LABEL org.opencontainers.image.url=https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite
LABEL org.opencontainers.image.source=https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite
LABEL org.opencontainers.image.title="odt"
LABEL org.opencontainers.image.description="Oligo Designer Toolsuite"
LABEL org.opencontainers.image.licenses=MIT
LABEL org.opencontainers.image.version=$ODT_REF

# --- Set up Python environment and install dependencies via conda ---
RUN --mount=source=environment.yaml,target=/tmp/env.yaml \
    --mount=type=cache,target=/opt/conda/pkgs,uid=$MAMBA_USER_ID \
    micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba install -y -n base -c conda-forge git

# --- Install ODT ---
# activate conda environment to use pip during build
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG ODT_REF
RUN --mount=type=cache,target=/home/$MAMBA_USER/.cache/pip,uid=$MAMBA_USER_ID \
    pip install "oligo-designer-toolsuite @ git+https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite.git@${ODT_REF}"
# ensure that docker/apptainer exec use the conda environment
# not recommended, but best solution: https://micromamba-docker.readthedocs.io/en/latest/advanced_usage.html#on-docker-exec
ENV PATH=/opt/conda/bin:$PATH

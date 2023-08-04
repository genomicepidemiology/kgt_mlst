FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY kmergenetyper.yml .

RUN conda env create -f kmergenetyper.yml

# Override default shell and use bash
SHELL ["conda", "run", "-n", "kmergenetyper", "/bin/bash", "-c"]

RUN kma -h

RUN pip install kmergenetyper

ENTRYPOINT ["conda", "run", "-n", "kmergenetyper"]

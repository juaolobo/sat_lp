FROM ubuntu:22.04

WORKDIR /latex-env

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Sao_Paulo

RUN apt update \
    && apt install -y texlive-base texlive-latex-base \
    texlive-plain-generic texlive-fonts-recommended \
    texlive-latex-recommended \
    texlive-latex-extra texlive-fonts-extra texlive-bibtex-extra \
    texlive-pictures texlive-science texlive-lang-portuguese \
    texlive-lang-english lmodern biber latexmk \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y make

ENTRYPOINT [ "/bin/bash" ]


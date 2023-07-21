
FROM --platform=linux/amd64 python:slim
RUN apt update && apt install procps -y
RUN apt install -y --no-install-recommends apt-utils
RUN apt -y install curl
RUN apt -y install libgomp1
LABEL MAINTAINER="J. Sebastian Paez - TalusBio"
COPY dist/*.whl /tmp/
WORKDIR /app
RUN WHEEL=$(ls /tmp/*.whl) && pip install ${WHEEL}

FROM python:3.11-slim

RUN apt-get update && \
    apt-get install -y infernal

RUN pip install -U pip && \
    pip install pybarrnap --no-cache-dir

CMD ["/bin/bash"]

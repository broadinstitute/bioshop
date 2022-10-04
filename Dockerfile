FROM python:bullseye
COPY . /app/src
RUN cd /app/src && \
    pip install . && \
    rm -rf /app

FROM python:bullseye
COPY . /app/src
RUN cd /app/src && \
    pip install . && \
    rm -rf /app
ENTRYPOINT ["/usr/local/bin/python", "/usr/local/bin/newt"]

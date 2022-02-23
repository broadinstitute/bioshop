FROM pytorch/pytorch:1.9.0-cuda10.2-cudnn7-runtime
COPY / /mostly/
RUN pip3 install -r /mostly/requirements.txt

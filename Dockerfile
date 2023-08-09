# golob/sepp:4.5.1


FROM --platform=amd64 debian:bullseye-slim

RUN export TZ=Etc/UTC
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get update && \
apt-get -y install tzdata && \
apt-get install -y \
    openjdk-11-jre \
    pigz \
    python3-pip \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ADD . /src/

RUN cd /src/ && \
python3 setup.py config && \
python3 setup.py install && \
mv /src/test /root/test/ && \
cd /root/ && rm -r /src/

WORKDIR /root/

ENTRYPOINT [ "/usr/local/bin/run_sepp.py" ] 
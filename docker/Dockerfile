# Galaxy - NRPS
#
# VERSION 0.1

FROM bgruening/galaxy-stable:18.01

MAINTAINER Loïc Couderc, loic.couderc@univ-lille1.fr

ENV GALAXY_CONFIG_BRAND nrps
ENV USER_PASSWORD nrps

COPY ./docker/job_conf.xml /etc/galaxy/job_conf.xml
COPY ./docker/update_local_slots_and_start_galaxy /usr/bin/update_local_slots_and_start_galaxy
COPY ./docker/welcome.html /etc/galaxy/web/welcome.html
#COPY ./galaxy/local_tools /local_tools

#Fix NRPSpredictor issue: https://github.com/roettig/NRPSpredictor2/issues/1
RUN sudo apt-get update && sudo apt-get install libc6-i386

# Mark folders as imported from the host.
VOLUME ["/export/", "/data/", "/var/lib/docker"]

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

# Add users
RUN service postgresql start && \
    . $GALAXY_VIRTUAL_ENV/bin/activate && \
    /galaxy_venv/bin/python /usr/local/bin/create_galaxy_user.py --user valerie.leclere@univ-lille1.fr --password $USER_PASSWORD --username valerie; \
    /galaxy_venv/bin/python /usr/local/bin/create_galaxy_user.py --user matthieu.duban@univ-lille1.fr --password $USER_PASSWORD --username matthieu;

# Add workflows
RUN mkdir -p $GALAXY_ROOT/workflows
COPY ./docker/workflows/* $GALAXY_ROOT/workflows/
RUN ls -la $GALAXY_ROOT/workflows/
RUN startup_lite && \
    sleep 60 && \
    . $GALAXY_VIRTUAL_ENV/bin/activate && \
    workflow-install --workflow_path $GALAXY_ROOT/workflows/ -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \
    workflow-install --workflow_path $GALAXY_ROOT/workflows/ -g http://localhost:8080 -u valerie.leclere@univ-lille1.fr -p $USER_PASSWORD && \
    workflow-install --workflow_path $GALAXY_ROOT/workflows/ -g http://localhost:8080 -u matthieu.duban@univ-lille1.fr -p $USER_PASSWORD

CMD ["/usr/bin/update_local_slots_and_start_galaxy"]

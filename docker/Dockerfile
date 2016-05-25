FROM ubuntu:latest

MAINTAINER Etienne Thevenot (etienne.thevenot@cea.fr)

# Install Ansible
RUN apt-get update
RUN apt-get -y install software-properties-common
RUN apt-add-repository ppa:ansible/ansible
RUN apt-get -y install ansible

# Run ansible role
ADD ansible /files/ansible
WORKDIR /files/ansible
RUN ansible-playbook multivariate.yml -c local

# Clean up
# TODO Remove also ansible stuff (ansible, Python, etc.)
RUN apt-get clean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/ /tmp/* /var/tmp/*

# Define Entry point script
ENTRYPOINT ["/files/multivariate/multivariate_wrapper.R"]

FROM ubuntu:14.04

MAINTAINER Etienne Thevenot (etienne.thevenot@cea.fr)

# Update package database
RUN apt-get update

# Install Ansible
RUN apt-get install -y git
RUN git clone git://github.com/ansible/ansible.git --recursive
WORKDIR ./ansible
RUN apt-get install -y python python-setuptools
RUN apt-get install -y gcc
RUN apt-get install -y python-dev
RUN apt-get install -y libgmp-dev
RUN easy_install pip
RUN pip install paramiko PyYAML Jinja2 httplib2 six
RUN apt-get install -y libffi-dev
RUN apt-get install -y libssl-dev
RUN pip install 'requests[security]'
RUN apt-get install -y make
RUN make all install
RUN ansible --version
RUN mkdir /etc/ansible
RUN echo '[local]\nlocalhost\n' > /etc/ansible/hosts

# Run ansible role
ADD ansible /files/ansible
WORKDIR /files/ansible
RUN ansible-playbook multivariate.yml -c local

# Clean up
# TODO Remove also ansible stuff (ansible, Python, etc.)
RUN apt-get clean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/ /tmp/* /var/tmp/*

# Define Entry point script
ENTRYPOINT ["/files/multivariate/multivariate_wrapper.R"]

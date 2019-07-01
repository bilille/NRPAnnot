## SYNOPSIS

Annotation of Nonribosomal Peptides in bacterial genomes with Galaxy.

This workflow uses the following tools for the identification and analysis of synthetase gene clusters:
* [antiSMASH3](https://dx.doi.org/10.1093%2Fnar%2Fgkv437)
* [NaPDoS](https://dx.doi.org/10.1371%2Fjournal.pone.0034064)

A comprehensive report is built with results from NaPDoS and the best ranked results of antiSMASH.

In practice, we provide two workflows :
* NRPSn: which takes a nucleic acid file as input (genbank/embl/fasta format) and performs the identification and analysis of synthetase gene clusters.
* NRPSp: which takes a protein sequence file as input (fasta format) directly with the proteins encoded by nrps genes from a cluster and perform the analysis of these synthetases.

## PREREQUISITES

You will need Docker installed on your execution environment.
See this [link](https://docs.docker.com/install/) for more details.

Then, add your user to docker group:
```bash
sudo usermod -aG docker your-user
## Logout and login again for this change to take effect
```

## INSTALL

```bash

## pull docker images (to avoid the loading time when needed by galaxy)
sudo docker pull antismash/standalone:3.0.5-2

## create export dir
sudo mkdir /export
sudo chown $USER /export

## clone the current repository
git clone https://github.com/bilille/NRPAnnot.git

## copy local_tools directory to $HOME/local_tools
cp -r NRPAnnot/galaxy/local_tools $HOME/

## run docker (restart automatically after reboot with --restart=always)
docker run -p 8000:80 --privileged=true -e DOCKER_PARENT=True -v /var/run/docker.sock:/var/run/docker.sock -v $HOME/local_tools:/local_tools -e "GALAXY_CONFIG_CONDA_AUTO_INSTALL=True" -v /export/:/export/ -e GALAXY_CONFIG_TOOL_CONFIG_FILE=config/tool_conf.xml.sample,config/shed_tool_conf.xml.sample,/local_tools/my_tools.xml -e "GALAXY_LOGGING=full" lcouderc/nrps-galaxy

## First run can takes some times due to the initialization of the container
## Open your favorite browser and open : http://127.0.0.1:8000/
## Log in with login `admin` and password `admin` to access the worklows.

```

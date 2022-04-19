# Introduction to containers (singularity)

**Table of contents**

<!-- vscode-markdown-toc -->
* 1. [Containers, images, ...  ?](#Containersimages...)
* 2. [What containers are not](#Whatcontainersarenot)
* 3. [Why do you (and don't) need containers](#Whydoyouanddontneedcontainers)
* 4. [References](#References)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

<HR>

**Aim of this session**

In this section we provide a very short introduction to the general idea behind containers. We briefly discuss the situations when containers can be a good choice for your problem and when other solutions might be more fitting. During this training event, we focus on one containerisation technology: Singularity, although we will cite other widespread container technologies such as Docker. You are most likely to see these two used in the research environment, with Singularity becoming the technology of choice for High Performance Computing development (although Docker is not going from that space anytime soon).


##  1. <a name='Containersimages...'></a>Containers, images, ...  ?

You may hear people say they are "running an image" or "running a container". These terms are often used interchangeably, although, they can mean different things - container being the running image.
The most important feature of containers, and where their real strength comes from, is that unlike "regular" applications, they can and often do perform all their work in isolation from their host OS.

Your containers do not have to know what the rest of your OS is up to. They don't even have to have access to the same files as your host OS, or share the same network  (again it is possible to achieve that). Containers put a layer between your existing host filesystem and whatever you are running inside them.

##  2. <a name='Whatcontainersarenot'></a>What containers are not

You will often hear the expression that "containers are like VMs", or "like VMs, but lighter". This may make sense on the surface. 

At the end of the day, containers make use of virtualisation, *BUT* a different kind of virtualisation. There are fewer moving components in the case of containers and the end result might be the same for the end user.

Containers remove a lot of components of virtual machines though: **they do not virtualise the hardware**, they do not have to contain a fully-fledged guest OS to operate. They have to **rely on the host OS** instead. VMs sit on top of the underlying hardware, whereas containers sit on top of the host OS.

![Containers vs VMs](https://gitlab.com/ska-telescope/src/ska-src-training-containers/-/raw/main/Talks/containers-intro/media/Virtual-Machine-and-Container-Deployments.png)
*Source [1](#ref1).*


For a more in-depth explanation of the differences between VMs and containers, please see this website by the [IBM Cloud Team](https://www.ibm.com/cloud/blog/containers-vs-vms).

##  3. <a name='Whydoyouanddontneedcontainers'></a>Why do you (and don't) need containers

**Yes, containers:**

- Containers will provide a reproducible work environment.
- They go beyond just sharing your code: you provide a fully-working software with all its required dependencies (modules, libraries, etc.).
- You can build self-contained images that meet the particular needs of your project. No need to install software "just in case", or install something to be used just once.
- You are no longer tied to the software and library versions installed on your host system.
- Need python3, but only python2 is available? There is an image for that.


**No, thanks:**

- Your software still depends on hardware you run it on - make sure your results are consistent across different hardware architectures.
- Not the best for sharing large amounts of data (same as you wouldn't use git to share a 10GB file).
- Additional safety concerns, as e.g. Docker gives extra power to the user "out of the box". There is potential to do some damage to the host OS by an inexperienced or malicious user if your containerisation technology of choice is not configured or used properly.




##  4. <a name='References'></a>References

1. <a name="ref1">NIST Special Publication 800-190, Application Container Security Guide - Scientific Figure on ResearchGate</a>. Available from: https://www.researchgate.net/figure/Virtual-Machine-and-Container-Deployments_fig1_329973333 [accessed 28 Jan, 2022]
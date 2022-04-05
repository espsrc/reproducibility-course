# Introduction to tools for repository management with `git`.

<HR>

**Table of contents**

<!-- vscode-markdown-toc -->
* 1. [Tools for `git`](#Toolsforgit)
* 2. [GIT desktop client](#GITdesktopclient)
* 3. [GIT from command line](#GITfromcommandline)
	* 3.1. [Git for Linux](#GitforLinux)
	* 3.2. [Git for MacOSX](#GitforMacOSX)
	* 3.3. [Git for Windows](#GitforWindows)
* 4. [Access to VM resources](#AccesstoVMresources)
	* 4.1. [How to access from a Web browser](#HowtoaccessfromaWebbrowser)
	* 4.2. [How to access from Remote Desktop for Windows](#HowtoaccessfromRemoteDesktopforWindows)
	* 4.3. [How to access from Remote Desktop for Linux](#HowtoaccessfromRemoteDesktopforLinux)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

<HR>

##  1. <a name='Toolsforgit'></a>Tools for `git`

As you know,  Git is a popular version control system that helps developers, writers, or anyone that requires keeping versions of their files to manage them and track changes. 

Many people use it mainly in the terminal or command line, although there are several applications that allow us to manage the use of our repository from completely visual tools, for example in the previous session we saw github from a website, but we also have the option of installing a git client for our computer as we will see later.

##  2. <a name='GITdesktopclient'></a>GIT desktop client

As we know, Git is necessary when it comes to managing collaborative development projects. Although, it also has a high learning curve. Thus, to make it easier for newcomers, developers, general public, etc.  have created Git Graphical User Interface clients for various platforms.

Here we recommend the standard github client called GitHub Desktop which works on Mac, Windows and Linux, is free and allows you to manage repo cloning, push-pull commands, merge conflicts, versioning, etc.

Here is the link to download: https://desktop.github.com/

Once downloaded, install it and you can start to use the application:

![first screen](./media/training1.1.png)

Then log-in from this application with GitHub

![2nd screen](./media/training1.2.png)

Finally you can start playing with the client of `git` from your personal computer desktop.

![3rd screen](./media/training1.3.png)

##  3. <a name='GITfromcommandline'></a>GIT from command line

There are many different ways to use Git. There are the original command-line tools, and there are many graphical user interfaces of different capabilities. The most widespread option for using Git is on the command line. 

For one thing, the command line is the only place where you can run all of Git's commands - most GUIs implement only a partial subset of Git's functionality for simplicity.

Also, while your choice of graphical client is a matter of personal taste, all users will have command-line tools installed and available on their system, especially if they use Linux or Mac (they already include it by default), and in the case of Windows, you can use some desktop client that adds command-line support.


###  3.1. <a name='GitforLinux'></a>Git for Linux

Depending on your Linux distribution, you can use the following options.

**Debian / Ubuntu (with apt-get)**

From your shell, install Git using apt-get:

```
sudo apt-get update
sudo apt-get install git
```

Verify the installation was successful by typing `git --version`

**Fedora (dnf or yum)**

From your shell, install Git using dnf (or yum, on older versions of Fedora and CentOS):

```
sudo dnf install git or sudo yum install git
```

Verify the installation was successful by typing `git --version`

###  3.2. <a name='GitforMacOSX'></a>Git for MacOSX

There are several ways to install Git on a Mac. If you’ve installed XCode (or it’s Command Line Tools), Git may already be installed. To find out, open a terminal and enter `git --version`. If the command returns the git version, then it is installed, otherwise you can install Xcode from the App Store or use the other methods below.

**Install Git on a Mac is via the stand-alone installer**

- Download the latest Git for Mac installer here and then install it.
- Open a terminal and verify the installation was successful by typing `git --version`.

**Install Git with Homebrew**

First you need to have HomeBrew installed, after that follow the next:

- Open your terminal and install Git using Homebrew: brew install git.
- Verify the installation was successful by typing `git --version`.

###  3.3. <a name='GitforWindows'></a>Git for Windows

Download the latest [Git for Windows](https://www.google.com/url?q=https://git-for-windows.github.io/&sa=D&source=editors&ust=1642764538497386&usg=AOvVaw0EAVw71DR5KUNL54vma1Zz) installer and then install it.

Then, open a Command Prompt (or Git Bash if during installation you elected not to use Git from the Windows Command Prompt), and verify the installation was successful by typing git --version and checking it’s showing a version greater or equal than 2.34.

##  4. <a name='AccesstoVMresources'></a>Access to VM resources

In order to make the most efficient use of this course we have created a Virtual Machine that has all the software resources that will be used in the reproducibility course. For this session we have enabled a workspace and a virtual desktop with git command line, in order to be able to do all the training of this session from a window (called terminal) that we have enabled.

To use these resources it is only necessary to have a web browser and to have received the access credentials to the Virtual Machine. 

###  4.1. <a name='HowtoaccessfromaWebbrowser'></a>How to access from a Web browser

To access open a tab in your browser and click on the following address: https://spsrc14.iaa.csic.es/guacamole/#/ 

![1 screen](./media/training2.1.png)

Use your username and password that you should have in the credentials and access email that has been sent for this second session. 

**:warning: :warning:  Once you have done this you will be asked again for your username and password, enter it by hand, for security reasons copy/paste has been disabled.**

![2 screen](./media/training2.2.png)

Once this is done you will see a linux desktop screen. 

![1 screen](./media/training2.3.png)

If you click on the 'Terminal' icon, you can follow the instructions for this session within the next steps.

![1 screen](./media/training2.4.png)


###  4.2. <a name='HowtoaccessfromRemoteDesktopforWindows'></a>How to access from Remote Desktop for Windows

This option is only if you are working from a Windows computer. To access the VM for the course, open the "Remote Desktop" (or "Escritorio Remoto") application.

Once you have done this, enter the following name where it says "Computer": 

> spsrc14.iaa.csic.es:18020

*This is the host and port to access with Remote Desktop*

![1 screen](./media/training3.2.png)

Click on connect and it will ask for your credentials. When you finish connecting (accept all the confirmations), the working desktop and the access terminal for `git` will appear.

![1 screen](./media/training2.3.png)

If you click on the 'Terminal' icon, you can follow the instructions for this session within the next steps.

![1 screen](./media/training2.4.png)


###  4.3. <a name='HowtoaccessfromRemoteDesktopforLinux'></a>How to access from Remote Desktop for Linux

This option is only if you are working from a Linux computer. To access the VM for the course, install REMMINA following this [tutorial](https://remmina.org/how-to-install-remmina/).

Once you have done this, enter the following name where it says "Server", change the protocol to "RDP" and add your credentials (user and password) to connect directly: 

![1 screen](./media/training4.1.png)

After that you will see:

![1 screen](./media/training2.3.png)

If you click on the 'Terminal' icon, you can follow the instructions for this session within the next steps.

![1 screen](./media/training2.4.png)